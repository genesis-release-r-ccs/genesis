!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_morph_mod
!> @brief   Moprhing simulation
!! @authors Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_morph_mod

  use at_energy_mod
  use at_energy_str_mod
  use at_energy_restraints_mod
  use at_dynvars_str_mod
  use at_dynvars_mod
  use at_ensemble_str_mod
  use at_morph_str_mod
  use at_restraints_str_mod
  use at_boundary_str_mod
  use at_pairlist_str_mod
  use at_pairlist_mod
  use at_enefunc_str_mod
  use at_enefunc_restraints_mod
  use at_output_str_mod
  use at_output_mod
  use molecules_str_mod
  use fileio_morph_mod
  use fileio_grotop_mod
  use fileio_control_mod
  use timers_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif

  implicit none
  private

  ! structures
  type, public :: s_morph_info
    integer          :: method              = MorphMinMethodSD
    integer          :: ncycles             =    100
    integer          :: iterations          =    100
    integer          :: iterations_2nd      =     10
    integer          :: eneout_period       =     10
    integer          :: crdout_period       =      0
    integer          :: velout_period       =      0
    integer          :: rstout_period       =      0
    integer          :: stoptr_period       =      0
    integer          :: nbupdate_period     =     10
    integer          :: morph_group         =      0
    integer          :: morph_bb_group      =      0
    integer          :: morph_sc_group      =      0
    real(wp)         :: morph_spring        = 5.0_wp
    real(wp)         :: morph_linear        = 1.0_wp
    real(wp)         :: morph_min_rmsd      = 0.5_wp
    real(wp)         :: morph_coef          = 1.5_wp
    real(wp)         :: morph_spring_max    = 20.0_wp
    logical          :: verbose             = .false.
  end type s_morph_info

  integer,  parameter :: metric = 10
  real(wp), parameter :: max_spring = 1.0E7_wp

  ! subroutines
  public  ::  show_ctrl_morph
  public  ::  read_ctrl_morph
  public  ::  setup_morph
  public  ::  setup_enefunc_morph_in
  public  ::  run_morph
  private ::  get_morph_in
  private ::  calcrmsd_morph
  private ::  morph_minimize_sd
  private ::  morph_minimize_lbfgs
  private ::  morph_minimize_steps
  private ::  check_morph

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_morph
  !> @brief        show DYNAMICS section usage
  !! @authors      CK
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "morph"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_morph(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('morph')

        write(MsgOut,'(A)') '[MORPH]'
        write(MsgOut,'(A)') 'method        = SD       # Method of Minimization'
        write(MsgOut,'(A)') 'ncycles       = 100       # number of Morphing cycles'
        write(MsgOut,'(A)') 'iterations    = 100       # number of iterations at each cycle'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') '# crdout_period = 0         # coordinates output period'
        write(MsgOut,'(A)') '# rstout_period = 0         # restart output period'
        write(MsgOut,'(A)') '# nbupdate_period = 10      # nonbond update period'
        write(MsgOut,'(A)') '# morph_spring = 5.0      # spring constant for morph'
        write(MsgOut,'(A)') '# morph_linear = 1.0      # linear constant for morph'
        write(MsgOut,'(A)') '# morph_min_rmsd = 1.5      # criterion minimum rmsd for morph'
        write(MsgOut,'(A)') '# morph_coef = 1.5      # coef for morph'
        write(MsgOut,'(A)') '# morph_spring_max = 20.0      # criterion spring constant for morph'
        write(MsgOut,'(A)') ' '


      end select

    else

      select case (run_mode)

      case ('morph')

        write(MsgOut,'(A)') '[MORPH]'
        write(MsgOut,'(A)') 'ncycles       = 100       # number of Morphing cycles'
        write(MsgOut,'(A)') 'iterations    = 100       # number of iterations at each cycle'
        write(MsgOut,'(A)') 'eneout_period = 10        # energy output period'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_morph

  !======1=========2=========3=========4=========5=========6=========7=========8 !
  !  Subroutine    read_ctrl_morph
  !> @brief        read DYNAMICS section in the control file
  !! @authors      TM, JJ
  !! @param[in]    handle     : unit number of control file
  !! @param[out]   morph_info : MORPH section in control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_morph(handle, morph_info)

    ! parameters
    character(*),            parameter     :: Section = 'Morph'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_morph_info),      intent(inout) :: morph_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'method',           &
                               morph_info%method, MorphMinMethodTypes)
    call read_ctrlfile_integer(handle, Section, 'ncycles',          &
                               morph_info%ncycles)
    call read_ctrlfile_integer(handle, Section, 'iterations',       &
                               morph_info%iterations)
    call read_ctrlfile_integer(handle, Section, 'iterations_2nd',   &
                               morph_info%iterations_2nd)
    call read_ctrlfile_integer(handle, Section, 'eneout_period',    &
                               morph_info%eneout_period)
    call read_ctrlfile_integer(handle, Section, 'crdout_period',    &
                               morph_info%crdout_period)
    call read_ctrlfile_integer(handle, Section, 'rstout_period',    &
                               morph_info%rstout_period)
    call read_ctrlfile_integer(handle, Section, 'nbupdate_period',  &
                               morph_info%nbupdate_period)
    call read_ctrlfile_integer(handle, Section, 'morph_group',      &
                               morph_info%morph_group)
    call read_ctrlfile_real   (handle, Section, 'morph_spring',     &
                               morph_info%morph_spring)
    call read_ctrlfile_real   (handle, Section, 'morph_spring_max', &
                               morph_info%morph_spring_max)
    call read_ctrlfile_real   (handle, Section, 'morph_linear',     &
                               morph_info%morph_linear)
    call read_ctrlfile_real   (handle, Section, 'morph_min_rmsd',   &
                               morph_info%morph_min_rmsd)
    call read_ctrlfile_real   (handle, Section, 'morph_coef',       &
                               morph_info%morph_coef)
    call read_ctrlfile_logical(handle, Section, 'verbose',          &
                               morph_info%verbose)
    call read_ctrlfile_integer(handle, Section, 'morph_sidechain',  &
                               morph_info%morph_sc_group)
    call read_ctrlfile_integer(handle, Section, 'morph_backbone',   &
                               morph_info%morph_bb_group)

    call end_ctrlfile_section(handle)

    ! error check
    !
    if (morph_info%crdout_period > 0) then
      if (mod(morph_info%ncycles, morph_info%crdout_period) /= 0) then
        call error_msg('Read_Ctrl_Morph> mod(ncycles, crdout_period) is not ZERO')
      end if
    end if

    if (morph_info%eneout_period > 0) then
      if (mod(morph_info%iterations, morph_info%eneout_period) /= 0) then
        call error_msg('Read_Ctrl_Morph> mod(iterations, eneout_period) is not ZERO')
      end if
    end if

    if (morph_info%rstout_period > 0) then
      if (mod(morph_info%ncycles, morph_info%rstout_period) /= 0) then
        call error_msg('Read_Ctrl_Morph> mod(ncycles, rstout_period) is not ZERO')
      end if
    end if

    if (morph_info%nbupdate_period > 0) then
      if (mod(morph_info%iterations, morph_info%nbupdate_period) /= 0) then
        call error_msg('Read_Ctrl_Morph> mod(iterations, nbupdate_period) is not ZERO')
      end if
    end if

    if (morph_info%morph_group <= 0) then
      call error_msg('Read_Ctrl_Morph> morph_group should be set.')
    end if

    if (morph_info%morph_spring_max <= morph_info%morph_spring) then
      call error_msg('Read_Ctrl_Morph> criterion of spring constant should be set.')
    end if

    if (morph_info%method == MorphMinMethodSTEPS) then 
      if (morph_info%morph_bb_group  == 0) &
         morph_info%morph_group = morph_info%morph_bb_group
      if (morph_info%morph_sc_group <= 0) then
        call error_msg('Read_Ctrl_Morph> morph_sc should be set.')
      end if
    end if


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Morph> Parameters of Morph simulation'

      write(MsgOut,'(A20,A10)')                                  &
            '  method          = ', MorphMinMethodTypes(morph_info%method)

      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  ncycles         = ', morph_info%ncycles,          &
            '  iterations      = ', morph_info%iterations
      if (morph_info%method == MorphMinMethodSTEPS) then 
        write(MsgOut,'(A20,I10)')                                &
            '  iterations_2nd  = ', morph_info%iterations_2nd
      endif
      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  eneout_period   = ', morph_info%eneout_period,    &
            '  rstout_period   = ', morph_info%rstout_period
      write(MsgOut,'(A20,I10,A20,I10)')                          &
            '  crdout_period   = ', morph_info%crdout_period,    &
            '  nbupdate_period = ', morph_info%nbupdate_period
      write(MsgOut,'(A20,I10)')                                  &
            '  morph_group     = ', morph_info%morph_group 

      if (morph_info%method == MorphMinMethodSTEPS) then 
        write(MsgOut,'(A20,I10,A20,I10)')                        &
            '  morph_backbone  = ', morph_info%morph_bb_group,   &
            '  morph_sidechain = ', morph_info%morph_sc_group 
      endif

      write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
            '  morph_spring    = ', morph_info%morph_spring,     &
            '  morph_linear    = ', morph_info%morph_linear
      write(MsgOut,'(A20,F10.3,A20,F10.3)')                      &
            '  morph_min_rmsd  = ', morph_info%morph_min_rmsd,   &
            '  morph_coef      = ', morph_info%morph_coef  

    end if

    return

  end subroutine read_ctrl_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  setup_morph
  !> @brief        setup morph
  !> @author       CK
  !! @param[in]    morph_info [type]
  !! @param[out]   morph [type]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_morph(morph_info, morph_in, morph, dynvars, molecule,  &
                         enefunc, restraints)

    ! formal arguments
    type(s_morph_info),       intent(in)    :: morph_info
    type(s_morph_in),         intent(in)    :: morph_in
    type(s_morph),            intent(inout) :: morph
    type(s_dynvars),  target, intent(in)    :: dynvars
    type(s_molecule), target, intent(in)    :: molecule
    type(s_enefunc),  target, intent(inout) :: enefunc
    type(s_restraints),       intent(in)    :: restraints

    integer                   :: i, j, nselect, nwork, n3
    logical                   :: flag
    real(wp)                  :: dx, dy, dz, dr
    real(wp)                  :: dx_t, dy_t, dz_t, dr_t
    real(wp)                  :: rmsd
    real(wp),         pointer :: coord(:,:),coord_ref(:,:)


    coord      => dynvars%coord
    coord_ref  => molecule%atom_refcoord
    
    call get_morph_in(morph_in, enefunc)

    ! setup morph
    !
    morph%method              = morph_info%method 
    morph%ncycles             = morph_info%ncycles 
    morph%iterations          = morph_info%iterations 
    morph%iterations_2nd      = morph_info%iterations_2nd
    morph%eneout_period       = morph_info%eneout_period
    morph%crdout_period       = morph_info%crdout_period
    morph%rstout_period       = morph_info%rstout_period
    morph%nbupdate_period     = morph_info%nbupdate_period
    morph%morph_min_rmsd      = morph_info%morph_min_rmsd
    morph%morph_coef          = morph_info%morph_coef  
    morph%verbose             = morph_info%verbose
    enefunc%morph_spring      = morph_info%morph_spring
    enefunc%morph_spring_max  = morph_info%morph_spring_max
    enefunc%morph_linear      = morph_info%morph_linear

    nselect = restraints%num_atoms(morph_info%morph_group)
    call alloc_morph(morph, MorphListRMSD, nselect)
    do i = 1, nselect
      morph%list_rmsd(i) = restraints%atomlist(i,morph_info%morph_group)
    end do

    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Morph> Setup Morph groups'
      write(MsgOut,'(A,I5,3A)') ' group = ',morph_info%morph_group
      write(MsgOut,'(A,I5)')   ' # of atoms = ', nselect
      write(MsgOut,'(A)')      ' atomlist: ' 
      do i = 1, nselect
        write(MsgOut,'(i7,$)') morph%list_rmsd(i)
        if (mod(i,10) == 0 .and. i /= nselect) then
          write(MsgOut,'(A)') ''
        end if
      end do
      write(MsgOut,'(A)') ''
    end if

    if (morph%method == MorphMinMethodSTEPS) then
      call alloc_morph(morph, MorphListBBSC, molecule%num_atoms)

      if (main_rank) then
        write(MsgOut,'(A)') 'Setup_Morph> Setup Morph Backbone groups'
        write(MsgOut,'(A,I5,3A)') ' group = ',morph_info%morph_bb_group
        write(MsgOut,'(A)') 'Setup_Morph> Setup Morph Sidechain groups'
        write(MsgOut,'(A,I5,3A)') ' group = ',morph_info%morph_sc_group
        write(MsgOut,'(A)') ''
      end if

      do i = 1, molecule%num_atoms
        morph%morph_atom_list(i) = MorphOther
        flag = .false.
        do j = 1, restraints%num_atoms(morph_info%morph_bb_group)
          if (i == restraints%atomlist(j,morph_info%morph_bb_group)) then
            flag = .true.
            exit
          end if
        end do

        if (flag) then

          morph%morph_atom_list(i) = MorphBackbone
        else

          do j = 1, restraints%num_atoms(morph_info%morph_sc_group)
            if (i == restraints%atomlist(j,morph_info%morph_sc_group)) then
              flag = .true.
              exit
            end if
          end do
          if (flag) morph%morph_atom_list(i) = MorphSidechain
        end if
      end do

    endif

    if (morph%method == MorphMinMethodLBFGS .or.  &
        morph%method == MorphMinMethodSTEPS ) then
      n3    = molecule%num_atoms*3
      nwork = 2*metric*n3 + 5*n3 + 11*metric*metric + 8*metric
      call alloc_morph(morph, MorphLBFGS, n3, nwork)
    end if

    call calcrmsd_morph(coord, coord_ref, molecule%mass, rmsd, morph)
    morph%morph_rmsd = rmsd
    morph%morph_rmsd_prev = rmsd
    morph%morph_drmsd =  &
       morph%morph_coef*morph%morph_rmsd/real(morph%ncycles, wp)
    if (main_rank) then
      write(MsgOut,'(A)') 'Setup_Morph> Initial Morphing check'
      write(MsgOut,'(A20,F10.3,A20,F10.3)')  &
            '     initial rmsd = ', morph%morph_rmsd,             &
            '  criterion drmsd = ', morph%morph_drmsd
      write(MsgOut,'(A)') ''
    end if

   return

  end subroutine setup_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  setup_enefunc_morph_in
  !> @brief        setup morph
  !> @author CK
  !! @param[in]    morph_in [type]
  !! @param[out]   enefunc [type]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_enefunc_morph_in(morph_in, enefunc)

    ! formal arguments
    type(s_morph_in),         intent(in) :: morph_in
    type(s_enefunc),  target, intent(inout) :: enefunc


    call get_morph_in(morph_in, enefunc)

   return

  end subroutine setup_enefunc_morph_in

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    calcrmsd_morph
  !> @brief        calculate rmsd
  !! @authors      CK
  !! @param[in]    coord: coordinate
  !! @param[in]    coord_ref: reference coordinate
  !! @param[in]    atomlist: atomlist
  !! @param[out]   rmsd
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine calcrmsd_morph(coord, coord_ref, mass, rmsd, morph)

    ! formal arguments
    real(wp),                intent(inout) :: coord(:,:)
    real(wp),                intent(in)    :: coord_ref(:,:)
    real(wp),                intent(in)    :: mass(:)
    real(wp),                intent(out)   :: rmsd
    type(s_morph),  target,  intent(inout) :: morph

    ! local variables
    real(wp)                 :: com_mov(1:3), com_ref(1:3)
    real(wp)                 :: tot_mass
    real(wp)                 :: dij(1:3)
    real(wp)                 :: sym_matrix(1:4,1:4)
    real(wp)                 :: rot_matrix(1:3,1:3), dadd(1:3)
    real(wp)                 :: eval(1:4), evec(1:4), lwork(1:11)
    real(wp)                 :: dref(1:3), dmov(1:3), dsub(1:3)
    integer                  :: nselect, iatm, i
    integer                  :: ierr

    real(wp), pointer        :: crd_mov(:,:), crd_rmov(:,:)
    real(wp), pointer        :: crd_rmov_fit(:,:), rmass(:)
    integer,  pointer        :: list(:)


    list         => morph%list_rmsd
    crd_mov      => morph%crd_mov
    crd_rmov     => morph%crd_rmov
    crd_rmov_fit => morph%crd_rmov_fit
    rmass        => morph%rmass

    nselect = size(list(:))
    com_mov(1:3) = 0.0_wp
    com_ref(1:3) = 0.0_wp
    tot_mass = 0.0_wp

    do i = 1, nselect
      iatm = list(i)
      com_mov(1:3)=com_mov(1:3)+mass(iatm)*coord(1:3,iatm)
      com_ref(1:3)=com_ref(1:3)+mass(iatm)*coord_ref(1:3,iatm)
      tot_mass = tot_mass + mass(iatm)
      rmass(i) = mass(iatm)
    end do
    com_mov(1:3) = com_mov(1:3)  / tot_mass
    com_ref(1:3) = com_ref(1:3) / tot_mass
    rmass(1:nselect) = rmass(1:nselect)/tot_mass

    ! symmetric matrix
    sym_matrix(1:4,1:4) = 0.0_wp
    do i = 1, nselect
      iatm = list(i)
      dref(1:3) = coord_ref(1:3,iatm) - com_ref(1:3)
      dmov(1:3) = coord(1:3,iatm)    - com_mov(1:3)
      dsub(1:3) = mass(iatm) * (dmov(1:3) - dref(1:3))
      dadd(1:3) = mass(iatm) * (dmov(1:3) + dref(1:3))
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

    ! eigenvalue and eigenvector
    !
#ifdef LAPACK
    call dsyev('V', 'U', 4, sym_matrix, 4, eval, lwork, 11, ierr)
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

    rmsd = 0.0_wp
    do i = 1, nselect
      iatm = list(i)
      crd_mov(1:3,i)      = coord(1:3,iatm)    - com_mov(1:3) 
      crd_rmov(1:3,i)     = coord_ref(1:3,iatm) - com_ref(1:3)
      crd_rmov_fit(1:3,i) = matmul(rot_matrix(1:3,1:3),crd_rmov(1:3,i))
      dij(1:3) = crd_mov(1:3,i) - crd_rmov_fit(1:3,i)
      rmsd = rmsd + rmass(i)*(dij(1)*dij(1) + dij(2)*dij(2) + dij(3)*dij(3))
    end do
    rmsd = sqrt(rmsd)

    return

  end subroutine calcrmsd_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine : run_morph
  !> @brief  entrance of the integrator
  !> @author CK
  !! @param[in] output:     information about output [str]
  !! @param[in] molecule:   information about molecules [str]
  !! @param[in] enefunc:    potential energy function [str]
  !! @param[in] dynvars:    dynamics varibales [str]
  !! @param[in] morph:      information about morph [str]
  !! @param[in] boundary:   boundary condition [str]
  !! @param[in] pairlist:   pairlist [str]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_morph(output, molecule, enefunc, dynvars, morph, &
                       boundary, pairlist)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_morph),            intent(inout) :: morph
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    integer                   :: i, ncycles, natom
    logical                   :: loopflg

    real(wp),         pointer :: coord(:,:), coord_ref(:,:)


    ! use pointers
    !
    ncycles    =  morph%ncycles
    natom      =  molecule%num_atoms
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref

    ! Open output files
    !
    call open_output(output)
    ! write 0 step
    call check_morph(dynvars, molecule, morph, enefunc, loopflg)
    do i = 1, natom
      coord_ref(1,i) = coord(1,i) 
      coord_ref(2,i) = coord(2,i) 
      coord_ref(3,i) = coord(3,i) 
    end do
    dynvars%step = 0
    call output_morph(output, molecule, morph, boundary, &
                      enefunc, dynvars)

    ! Main loop of the Steepest descent method
    !
    i = 0
    loopflg = .true.
    do while(loopflg)
      i = i+1
      dynvars%step = i
      dynvars%iterations = 0

      select case (morph%method)

        case(MorphMinMethodSD)
       
          call morph_minimize_sd(output, molecule, enefunc, dynvars, morph, &
                                boundary, pairlist)
       
        case(MorphMinMethodLBFGS)
          call morph_minimize_lbfgs(output, molecule, enefunc, dynvars, morph, &
                                boundary, pairlist)

        case(MorphMinMethodSTEPS)
          call morph_minimize_steps(output, molecule, enefunc, dynvars, morph, &
                                boundary, pairlist)

      end select

      ! output trajectory
      !
      !
      call output_morph(output, molecule, morph, boundary, &
                                 enefunc, dynvars)

      call check_morph(dynvars, molecule, morph, enefunc, loopflg)
      if (i == ncycles) loopflg = .false.

    end do

    ! Close output files
    !
    call close_output(output)

    return

  end subroutine run_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  morph_minimize_sd
  !> @brief        entrance of the integrator
  !> @author       CK
  !! @param[in] output:   information about output [str]
  !! @param[in] molecule: information about molecules [str]
  !! @param[in] enefunc:  potential energy function [str]
  !! @param[in] dynvars:  dynamics varibales [str]
  !! @param[in] moorph:   information about moorph [str]
  !! @param[in] boundary: boundary condition [str]
  !! @param[in] pairlist: pairlist [str]
  !! @date   CK: 2012/11/13
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine morph_minimize_sd(output, molecule, enefunc, dynvars, morph, &
                               boundary, pairlist)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_morph),    target, intent(inout) :: morph
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)                  :: energy_ref, delta_energy
    real(wp)                  :: delta_rini, delta_rmax
    real(wp)                  :: rmsg, maxg, absg
    integer                   :: i, j, natom, nsteps
    real(wp),         pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),         pointer :: force_omp(:,:,:)
    real(wp),         pointer :: delta_r


    ! use pointers
    !
    natom      =  molecule%num_atoms
    coord      => dynvars%coord
    coord_ref  => dynvars%coord_ref
    force      => dynvars%force
    force_omp  => dynvars%force_omp
    delta_r    => morph%delta_r
    nsteps     =  morph%iterations

    delta_rmax = 0.0001_wp
    delta_rini = 0.00005_wp
    delta_r = delta_rini


    ! Compute energy of the initial structure
    !
    call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            .false.,                                  &
                            coord, dynvars%trans, dynvars%coord_pbc,  &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

    rmsg = 0.0_wp
    maxg = 0.0_wp
    do j = 1, natom
      rmsg = rmsg + force(1,j)**2 + force(2,j)**2 + force(3,j)**2
      absg = max(abs(force(1,j)),abs(force(2,j)),abs(force(3,j)))
      if(absg > maxg) maxg = absg
    end do
    rmsg = sqrt(rmsg/real(3*natom,wp))
    dynvars%rms_gradient = rmsg
    dynvars%max_gradient = maxg

    coord(1:3,1:natom) = coord(1:3,1:natom) + delta_r*force(1:3,1:natom)/rmsg

    ! Main loop of the Steepest descent method
    !
    do i = 1, nsteps

      dynvars%iterations = i
      ! save old energy and coordinates
      !
      coord_ref(1:3,1:natom) = coord(1:3,1:natom)
      energy_ref = dynvars%energy%total


      ! Compute energy and forces
      !
      call compute_energy(molecule, enefunc, pairlist, boundary,      &
                          mod(i, morph%eneout_period) == 0,           &
                          .false.,                                    &
                          coord, dynvars%trans, dynvars%coord_pbc,    &
                          dynvars%energy, dynvars%temporary,          &
                          force, force_omp, dynvars%virial,           &
                          dynvars%virial_extern)

      delta_energy = dynvars%energy%total - energy_ref
      if (delta_energy > 0) then
        delta_r = 0.5_wp*delta_r
      else if (abs(delta_energy) < morph%minimize_cutoff) then
        !
        ! energy conversed
        !
        return
      else
        delta_r = min(delta_rmax, 1.2_wp*delta_r)
      end if

      rmsg = 0.0_wp
      maxg = 0.0_wp
      do j = 1, natom
        rmsg = rmsg + force(1,j)**2 + force(2,j)**2 + force(3,j)**2
        absg = max(abs(force(1,j)),abs(force(2,j)),abs(force(3,j)))
        if(absg > maxg) maxg = absg
      end do
      rmsg = sqrt(rmsg/real(3*natom,wp))
      dynvars%rms_gradient = rmsg
      dynvars%max_gradient = maxg

      coord(1:3,1:natom) = coord(1:3,1:natom) + delta_r*force(1:3,1:natom)/rmsg

      if (mod(i,morph%eneout_period) == 0) then
        dynvars%total_pene = dynvars%energy%total
        call output_dynvars(output, enefunc, dynvars)
      endif

      ! update nonbond pairlist
      !
      if (morph%nbupdate_period > 0) then
        if (mod(i,morph%nbupdate_period) == 0 .and. real_calc) then

            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 dynvars%coord_pbc, pairlist)

        end if

      end if

    end do

    return

  end subroutine morph_minimize_sd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  morph_minimize_lbfgs
  !> @brief        entrance of the integrator
  !> @author       CK
  !! @param[in]    output:   information about output [str]
  !! @param[in]    molecule: information about molecules [str]
  !! @param[in]    enefunc:  potential energy function [str]
  !! @param[in]    dynvars:  dynamics varibales [str]
  !! @param[in]    moorph:   information about moorph [str]
  !! @param[in]    boundary: boundary condition [str]
  !! @param[in]    pairlist: pairlist [str]
  !! @date   CK: 2012/11/13
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine morph_minimize_lbfgs(output, molecule, enefunc, dynvars, morph, &
                                  boundary, pairlist)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_morph),    target, intent(inout) :: morph
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)               :: dsave(29)
    real(wp)               :: rmsg, maxg, absg
    real(wp)               :: energy_ref
    integer                :: i, j, ii, natom, nsteps
    integer                :: itime
    integer                :: isave(44)
    integer                :: iprint = -1
    character(len=60)      :: task, csave
    logical                :: lsave(4)

    real(wp),      pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),      pointer :: force_omp(:,:,:)
    real(wp),      pointer :: lower(:), upper(:)
    real(wp),      pointer :: vec(:), gradient(:), work_lbfgs(:)
    integer,       pointer :: list_bound(:), iwork_lbfgs(:)

    real(wp),    parameter :: factr  = 1.0d7, pgtol  = 1.0d-5

!    if (main_rank) then
!      iprint = 1
!    endif


    ! use pointers
    !
    natom       =  molecule%num_atoms
    coord       => dynvars%coord
    coord_ref   => dynvars%coord_ref
    force       => dynvars%force
    force_omp   => dynvars%force_omp
    vec         => morph%vec
    gradient    => morph%gradient
    work_lbfgs  => morph%work_lbfgs
    lower       => morph%lower
    upper       => morph%upper
    iwork_lbfgs => morph%iwork_lbfgs
    list_bound  => morph%list_bound

    nsteps     =  morph%iterations

    ii = 0
    do i = 1, natom
      ii = ii+1
      list_bound(ii) = 0
      vec(ii) = coord(1,i)
      ii = ii+1
      list_bound(ii) = 0
      vec(ii) = coord(2,i)
      ii = ii+1
      list_bound(ii) = 0
      vec(ii) = coord(3,i)
    end do

    task = 'START'

    ! Compute energy of the initial structure
    !
    itime = 0
    do while((task(1:2) .eq. 'FG' .or. task .eq. 'NEW_X' .or. &
             task .eq. 'START') .and. itime <= nsteps) 
      itime = itime + 1
      dynvars%iterations = itime

      call setulb(3*natom, metric, vec, lower ,upper, list_bound, energy_ref, &
                   gradient, factr, pgtol, work_lbfgs, iwork_lbfgs,           &
                   task, iprint, csave, lsave, isave, dsave)

      if (task(1:2) .eq. 'FG') then
        ii = 0
        do i = 1, natom
         ii = ii+1
         coord(1,i) = vec(ii) 
         ii = ii+1
         coord(2,i) = vec(ii) 
         ii = ii+1
         coord(3,i) = vec(ii) 
        end do
        coord_ref(1:3,1:natom) = coord(1:3,1:natom)
        call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                            .false.,                                  &
                            coord, dynvars%trans, dynvars%coord_pbc,  &
                            dynvars%energy, dynvars%temporary,        &
                            force, force_omp, dynvars%virial,         &
                            dynvars%virial_extern)

        energy_ref = dynvars%energy%total
        ii = 0
        do i = 1, natom
         ii = ii+1
         gradient(ii) = -force(1,i)
         ii = ii+1
         gradient(ii) = -force(2,i)
         ii = ii+1
         gradient(ii) = -force(3,i)
        end do
      end if

      rmsg = 0.0_wp
      maxg = 0.0_wp
      do j = 1, natom
        rmsg = rmsg + force(1,j)**2 + force(2,j)**2 + force(3,j)**2
        absg = max(abs(force(1,j)),abs(force(2,j)),abs(force(3,j)))
        if(absg > maxg) maxg = absg
      end do
      rmsg = sqrt(rmsg/real(3*natom,wp))
      dynvars%rms_gradient = rmsg
      dynvars%max_gradient = maxg

      if (mod(itime,morph%eneout_period) == 0) then
        dynvars%total_pene = dynvars%energy%total
        call output_dynvars(output, enefunc, dynvars)
      endif

      ! update nonbond pairlist
      !
      if (morph%nbupdate_period > 0) then
        if (mod(itime,morph%nbupdate_period) == 0 .and. real_calc) then

            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 dynvars%coord_pbc, pairlist)

        end if
      end if

    end do

    return

  end subroutine morph_minimize_lbfgs

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  morph_minimize_steps
  !> @brief        entrance of the integrator
  !> @author       CK
  !! @param[in]    output:   information about output [str]
  !! @param[in]    molecule: information about molecules [str]
  !! @param[in]    enefunc:  potential energy function [str]
  !! @param[in]    dynvars:  dynamics varibales [str]
  !! @param[in]    moorph:   information about moorph [str]
  !! @param[in]    boundary: boundary condition [str]
  !! @param[in]    pairlist: pairlist [str]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine morph_minimize_steps(output, molecule, enefunc, dynvars, morph, &
                                  boundary, pairlist)

    ! formal arguments
    type(s_output),           intent(inout) :: output
    type(s_molecule), target, intent(inout) :: molecule
    type(s_enefunc),          intent(inout) :: enefunc
    type(s_dynvars),  target, intent(inout) :: dynvars
    type(s_morph),    target, intent(inout) :: morph
    type(s_pairlist),         intent(inout) :: pairlist
    type(s_boundary),         intent(inout) :: boundary

    ! local variables
    real(wp)               :: energy_ref
    real(wp)               :: dsave(29)
    real(wp)               :: rmsg, maxg, absg
    real(wp)               :: spring_save
    integer                :: i, j, ii, natom
    integer                :: itime, istep
    integer                :: nstep_morph(2)
    integer                :: isave(44)
    integer                :: iprint = -1
    character(len=60)      :: task, csave
    logical                :: lsave(4)

    real(wp),      pointer :: coord(:,:), coord_ref(:,:), force(:,:)
    real(wp),      pointer :: force_omp(:,:,:)
    real(wp),      pointer :: lower(:), upper(:)
    real(wp),      pointer :: vec(:), gradient(:), work_lbfgs(:)
    integer,       pointer :: list_bound(:), iwork_lbfgs(:)

    real(wp),    parameter :: factr  = 1.0d7, pgtol  = 1.0d-5
    integer,     parameter :: move_steps(2) = (/-1, MorphBackbone/)


    ! use pointers
    !
    natom       =  molecule%num_atoms
    coord       => dynvars%coord
    coord_ref   => dynvars%coord_ref
    force       => dynvars%force
    force_omp   => dynvars%force_omp
    vec         => morph%vec
    gradient    => morph%gradient
    work_lbfgs  => morph%work_lbfgs
    lower       => morph%lower
    upper       => morph%upper
    iwork_lbfgs => morph%iwork_lbfgs
    list_bound  => morph%list_bound

    nstep_morph(1) = morph%iterations
    nstep_morph(2) = morph%iterations_2nd

    list_bound(1:3*natom) = 0
    istep = 0
    enefunc%morph_ene_flag = MorphEneFlagBB

    ! Compute energy of the initial structure
    !
    do while (istep < 2)
      ii = 0
      itime = 0
      istep = istep + 1
      do i = 1, natom
        ii = ii+1
        vec(ii) = coord(1,i)
        ii = ii+1
        vec(ii) = coord(2,i)
        ii = ii+1
        vec(ii) = coord(3,i)
      end do
      upper(1:3*natom) = vec(1:3*natom)
      lower(1:3*natom) = vec(1:3*natom)
      ii = 0
      do i = 1, natom
       if (morph%morph_atom_list(i) == move_steps(istep)) then
         ii = ii+1
         list_bound(ii) = 2
         ii = ii+1
         list_bound(ii) = 2
         ii = ii+1
         list_bound(ii) = 2
       end if
      end do

      if (istep == 2) then
        spring_save = enefunc%morph_spring
        enefunc%morph_spring = morph%morph_step_ratio
        enefunc%morph_ene_flag = MorphEneFlagSC
      end if

      task = 'START'
      do while((task(1:2) .eq. 'FG' .or. task .eq. 'NEW_X' .or. &
               task .eq. 'START') .and. itime <= nstep_morph(istep)) 
        itime = itime + 1
        dynvars%iterations = itime
     
        call setulb(3*natom, metric, vec, lower ,upper,                 &
                     list_bound, energy_ref,                            &
                     gradient, factr, pgtol, work_lbfgs, iwork_lbfgs,   &
                     task, iprint, csave, lsave, isave, dsave)
     
        if (task(1:2) .eq. 'FG') then
          ii = 0
          do i = 1, natom
            ii = ii+1
            coord(1,i) = vec(ii) 
            ii = ii+1
            coord(2,i) = vec(ii) 
            ii = ii+1
            coord(3,i) = vec(ii) 
          end do
          coord_ref(1:3,1:natom) = coord(1:3,1:natom)
          call compute_energy(molecule, enefunc, pairlist, boundary, .true., &
                              .false.,                                  &
                              coord, dynvars%trans, dynvars%coord_pbc,  &
                              dynvars%energy, dynvars%temporary,        &
                              force, force_omp, dynvars%virial,         &
                              dynvars%virial_extern)
     
          energy_ref = dynvars%energy%total
          dynvars%rms_gradient = rmsg
          ii = 0
          do i = 1, natom
            ii = ii+1
            gradient(ii) = -force(1,i)
            ii = ii+1
            gradient(ii) = -force(2,i)
            ii = ii+1
            gradient(ii) = -force(3,i)
          end do
        end if
     
        if (mod(itime,morph%eneout_period) == 0 .and. istep == 1) then
          dynvars%total_pene = dynvars%energy%total
          call output_dynvars(output, enefunc, dynvars)
        endif
     
        ! update nonbond pairlist
        !
        if (morph%nbupdate_period > 0) then
          if (mod(itime,morph%nbupdate_period) == 0 .and. real_calc) then
     
            call update_pairlist(enefunc, boundary, coord, dynvars%trans, &
                                 dynvars%coord_pbc, pairlist)
     
          end if
        end if
     
      end do
    end do
    enefunc%morph_spring = spring_save

    return

  end subroutine morph_minimize_steps

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  check_morph
  !> @brief        check_morph
  !> @author       CK
  !! @param[out]   dynvars [type]
  !! @param[out]   molecule [type]
  !! @param[out]   morph [type]
  !! @param[out]   enefunc [type]
  !> @date   CK: 2012/11/13
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_morph(dynvars, molecule, morph, enefunc, loopflg)

    ! formal arguments
    type(s_morph),           intent(inout) :: morph
    type(s_dynvars), target, intent(inout) :: dynvars
    type(s_enefunc), target, intent(inout) :: enefunc
    type(s_molecule),target, intent(in)    :: molecule
    logical,                 intent(inout) :: loopflg

    ! local variables
    real(wp)                 :: drmsd,rmsd
    real(wp),        pointer :: coord(:,:),coord_ref(:,:)


    coord      => dynvars%coord
    coord_ref  => molecule%atom_refcoord

    call calcrmsd_morph(coord, coord_ref, molecule%mass, rmsd, morph)

    morph%morph_rmsd = rmsd
    drmsd = morph%morph_rmsd_prev-morph%morph_rmsd
    morph%morph_rmsd_prev = morph%morph_rmsd

    if (main_rank) then
      write(MsgOut,'(A)') 'Check_morph> Morphing check'
      write(MsgOut,'(A16,F11.3,A16,F11.3,A16,F11.3)')  &
            '         rmsd = ', morph%morph_rmsd,             &
            '        drmsd = ', drmsd,                        &
            ' morph_spring = ', enefunc%morph_spring
    endif
    if (morph%morph_rmsd < morph%morph_min_rmsd) then
      loopflg = .false.
      return
    end if

    if (drmsd > morph%morph_drmsd) then
      enefunc%morph_spring = enefunc%morph_spring*morph%decrease
    else
      enefunc%morph_spring = enefunc%morph_spring*morph%increase
    end if
    if (enefunc%morph_spring > enefunc%morph_spring_max) then
      loopflg = .false.
      if (main_rank) then
        write(MsgOut,'(A)') 'Check_morph> morph_spring is greater than morph_spring_max'
        write(MsgOut,'(A16,F11.3,A16,F11.3,A16,F11.3)')  &
              '         rmsd = ', morph%morph_rmsd,             &
              '        drmsd = ', drmsd,                        &
              ' morph_spring = ', enefunc%morph_spring
      end if
      return
!      call error_msg('Check_Morph> morph_spring is greater than morph_spring_max')
    end if
    if (main_rank) then
      write(MsgOut,'(A22,F11.3)')  &
            ' after morph_spring = ', enefunc%morph_spring
      write(MsgOut,*) 
    end if
    if (enefunc%morph_spring > max_spring) then
      if (main_rank) then
        write(MsgOut,'(A)') 'Spring Exceeds Maximum: Stop'
      end if
      loopflg = .false.
    end if

    return

  end subroutine check_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine :  get_morph_in
  !> @brief        setup morph
  !> @author       CK
  !! @param[in]    morph_in [type]
  !! @param[out]   enefunc [type]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine get_morph_in(morph_in, enefunc)

    ! formal arguments
    type(s_morph_in),        intent(in)    :: morph_in
    type(s_enefunc), target, intent(inout) :: enefunc

    ! local variables
    integer                   :: i, nmorph_bb, nmorph_sc
    integer                   :: istart, iend


    if (enefunc%num_morph_bb + enefunc%num_morph_sc > 0) return

    enefunc%num_morph_bb = morph_in%num_morph_bb
    enefunc%num_morph_sc = morph_in%num_morph_sc

    call alloc_enefunc(enefunc, EneFuncMorph, enefunc%num_morph_bb,  &
                      enefunc%num_morph_sc)
    if (main_rank) then
      write(MsgOut,'(A)') 'Get_Morph_in> Morphing check'
      write(MsgOut,'(A20,I10,A20,I10)')  &
            '  morphing_bb     = ', enefunc%num_morph_bb, &
            '  morphing_sc     = ', enefunc%num_morph_sc
      write(MsgOut,'(A)') ' '
    end if

    nmorph_bb = 0
    do i = 1, morph_in%num_morph_bb
      nmorph_bb = nmorph_bb + 1
      enefunc%morph_list_bb(1,nmorph_bb) = &
             morph_in%morph_bb(nmorph_bb)%atom_idx1 
      enefunc%morph_list_bb(2,nmorph_bb) = &
             morph_in%morph_bb(nmorph_bb)%atom_idx2
      enefunc%morph_dist_bb(nmorph_bb)   = morph_in%morph_bb(nmorph_bb)%rmin
      enefunc%morph_dist_bb_other(nmorph_bb)   = morph_in%morph_bb(nmorph_bb)%rmin_other
    end do 
    nmorph_sc = 0
    do i = 1, morph_in%num_morph_sc
      nmorph_sc = nmorph_sc + 1
      enefunc%morph_list_sc(1,nmorph_sc) = &
             morph_in%morph_sc(nmorph_sc)%atom_idx1 
      enefunc%morph_list_sc(2,nmorph_sc) = &
             morph_in%morph_sc(nmorph_sc)%atom_idx2
      enefunc%morph_dist_sc(nmorph_sc)   = morph_in%morph_sc(nmorph_sc)%rmin
      enefunc%morph_dist_sc_other(nmorph_sc)   = morph_in%morph_sc(nmorph_sc)%rmin_other
    end do 

    call get_loop_index(enefunc%num_morph_bb, istart, iend)
    enefunc%istart_morph_bb = istart
    enefunc%iend_morph_bb   = iend

    call get_loop_index(enefunc%num_morph_sc, istart, iend)
    enefunc%istart_morph_sc = istart
    enefunc%iend_morph_sc   = iend
      
    return

  end subroutine get_morph_in

end module at_morph_mod
