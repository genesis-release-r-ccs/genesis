!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fitting_mod
!> @brief   module for fitting
!! @authors Yuji Sugita (YS), Norio Takase (NT), Chigusa Kobayashi (CK)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fitting_mod

  use fitting_str_mod
  use select_mod
  use select_atoms_mod
  use molecules_str_mod
  use messages_mod
  use fileio_control_mod
  use mpi_parallel_mod
  use constants_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_fit_info
    real(wp)   :: grid_size         = 1.0_wp
    integer    :: fitting_method    = FittingMethodNO
    integer    :: fitting_atom      = 1
    integer    :: num_grids         = 10
    logical    :: mass_weight       = .false.
    logical    :: force_no_fitting  = .false.
  end type s_fit_info

  public  :: show_ctrl_fitting
  public  :: read_ctrl_fitting
  public  :: read_ctrl_fitting_md
  public  :: setup_fitting
  public  :: run_fitting
  public  :: out_rmsd
  public  :: out_trrot
  public  :: fit_trrot
  public  :: fit_trans
  public  :: fit_zrot
  public  :: fit_xytr_zrot2
  public  :: transform

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_fitting
  !> @brief        show control parameters in FITTING section
  !! @authors      NT
  !! @param [in]   tool_name : tool name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_fitting(tool_name)

    ! formal arguments
    character(len=*), optional, intent(in) :: tool_name


    if (main_rank) then
      if (.not.present(tool_name)) then
        write(MsgOut,'(A)') '[FITTING]'
        write(MsgOut,'(A)') 'fitting_method = NO              # NO/TR+ROT/TR/TR+ZROT/XYTR/XYTR+ZROT'
        write(MsgOut,'(A)') 'fitting_atom   = 1               # atom group'
        write(MsgOut,'(A)') 'zrot_ngrid     = 10              # number of z-rot grids'
        write(MsgOut,'(A)') 'zrot_grid_size = 1.0             # z-rot grid size'
        write(MsgOut,'(A)') 'mass_weight    = NO              # mass-weight is not applied'

      else
        select case (tool_name)

        case ('qg') 
          write(MsgOut,'(A)') '# [FITTING]'
          write(MsgOut,'(A)') '# fitting_method = TR+ROT          # NO/TR+ROT/TR/TR+ZROT/XYTR/XYTR+ZROT'
          write(MsgOut,'(A)') '# fitting_atom   = 5               # atom group'
          write(MsgOut,'(A)') '# zrot_ngrid     = 10              # number of z-rot (ZROT) grids'
          write(MsgOut,'(A)') '# zrot_grid_size = 1.0             # z-rot (ZROT) grid size'
          write(MsgOut,'(A)') '# mass_weight    = NO              # mass-weight is not applied'

        case default
          write(MsgOut,'(A)') '[FITTING]'
          write(MsgOut,'(A)') 'fitting_method = NO              # NO/TR+ROT/TR/TR+ZROT/XYTR/XYTR+ZROT'
          write(MsgOut,'(A)') 'fitting_atom   = 1               # atom group'
          write(MsgOut,'(A)') 'zrot_ngrid     = 10              # number of z-rot grids'
          write(MsgOut,'(A)') 'zrot_grid_size = 1.0             # z-rot grid size'
          write(MsgOut,'(A)') 'mass_weight    = NO              # mass-weight is not applied'

        end select
      end if
      write(MsgOut,'(A)') ' '
    endif

    return

  end subroutine show_ctrl_fitting

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_fitting
  !> @brief        read FITTING section in the control file
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   fit_info : FITTING section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_fitting(handle, fit_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Fitting'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_fit_info),        intent(inout) :: fit_info


    ! read parameters 
    !

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'fitting_method', &
                               fit_info%fitting_method, FittingMethodTypes)

    call read_ctrlfile_integer(handle, Section, 'fitting_atom', &
                               fit_info%fitting_atom)

    call read_ctrlfile_integer(handle, Section, 'zrot_ngrid', &
                               fit_info%num_grids)

    call read_ctrlfile_real   (handle, Section, 'zrot_grid_size', &
                               fit_info%grid_size)
    call read_ctrlfile_logical(handle, Section, 'mass_weight',    &
                               fit_info%mass_weight)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Fitting> Parameters of Fitting'
      write(MsgOut,'(A20,A10)') &
                      '  fitting method  = ', &
                      FittingMethodTypes(fit_info%fitting_method)
      write(MsgOut,'(A20,A5,I0)') &
                      '  fitting atom    = ', 'group', fit_info%fitting_atom
      write(MsgOut,'(A20,I10)') &
                      '  z-rot # of grid = ', fit_info%num_grids
      write(MsgOut,'(A20,F10.3)') &
                      '  z-rot grid size = ', fit_info%grid_size
      if (fit_info%mass_weight) then
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =        yes'
      else
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =         no'
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_fitting

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_fitting_md
  !> @brief        read FITTING section in the control file
  !! @authors      NT, CK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   fit_info : FITTING section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_fitting_md(handle, fit_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Fitting'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_fit_info),        intent(inout) :: fit_info


    ! read parameters 
    !
    fit_info%fitting_method = FittingMethodTR_ROT

    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_type   (handle, Section, 'fitting_method', &
                               fit_info%fitting_method, FittingMethodTypes)

    call read_ctrlfile_integer(handle, Section, 'fitting_atom', &
                               fit_info%fitting_atom)

    call read_ctrlfile_logical(handle, Section, 'mass_weight',    &
                               fit_info%mass_weight)

    call read_ctrlfile_logical(handle, Section, 'force_no_fitting',    &
                               fit_info%force_no_fitting)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Fitting_MD> Parameters of Fitting'
      write(MsgOut,'(A20,A10)') &
                      '  fitting method  = ', &
                      FittingMethodTypes(fit_info%fitting_method)
      write(MsgOut,'(A20,A5,I0)') &
                      '  fitting atom    = ', 'group', fit_info%fitting_atom
      if (fit_info%mass_weight) then
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =        yes'
      else
        write(MsgOut,'(A30)')  &
                      '  mass_weight     =         no'
      end if
      if (fit_info%force_no_fitting) then
        write(MsgOut,'(A30)')  &
                      ' force_no_fitting =        yes'
      end if
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_fitting_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_fitting
  !> @brief        setup fitting information
  !! @authors      NT
  !! @param[in]    fit_info : FITTING section control parameters information
  !! @param[in]    sel_info : SELECTION section control parameters information
  !! @param[in]    molecule : molecule information
  !! @param[inout] fitting  : fitting information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine setup_fitting(fit_info, sel_info,  molecule, fitting)
  
    ! formal argments
    type(s_fit_info),        intent(in)    :: fit_info
    type(s_sel_info),        intent(in)    :: sel_info
    type(s_molecule),        intent(in)    :: molecule
    type(s_fitting),         intent(inout) :: fitting


    ! fitting method
    fitting%fitting_method = fit_info%fitting_method

    ! fitting atom
    if (fitting%fitting_method /= FittingMethodNo) then

      if (fit_info%fitting_atom > size(sel_info%groups)) &
        call error_msg('Setup_Fitting> fitting atom seleciton is out of range.')

      call select_atom(molecule, &
                       sel_info%groups(fit_info%fitting_atom), &
                       fitting%fitting_atom)

      write(MsgOut,'(A,I8)') 'Setup_Fitting> fitting atom count: ', &
           size(fitting%fitting_atom%idx)
      write(MsgOut,'(A)') ' '

    end if

    ! num grid
    fitting%num_grids = fit_info%num_grids

    ! grid size
    fitting%grid_size = fit_info%grid_size

    ! z-translation 
    fitting%z_translation = fitting%fitting_method < FittingMethodXYTR

    ! mass-wegiht
    fitting%mass_weight = fit_info%mass_weight

    return

  end subroutine setup_fitting

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    run_fitting
  !> @brief        run fitting
  !! @authors      YS, NT, CK
  !! @param[inout] fitting    : fitting information
  !! @param[in]    coord_ref  : reference coordinate data
  !! @param[in]    coord_mov  : moving coordinate data
  !! @param[out]   coord_fit  : fitted coordinate data
  !! @param[in]    mass_ref   : reference mass data (optional)
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_fitting(fitting,   &
                         coord_ref, &
                         coord_mov, &
                         coord_fit, &
                         mass_ref)

    ! formal arguments
    type(s_fitting),         intent(inout) :: fitting
    real(wp),                intent(in)    :: coord_ref(:,:)
    real(wp),                intent(in)    :: coord_mov(:,:)
    real(wp),                intent(out)   :: coord_fit(:,:)
    real(wp),      optional, intent(in)    :: mass_ref(:)

    ! local variables
    integer                  :: num_fit_atom, num_mov_atom
    real(wp),    allocatable :: mass(:)


    if (fitting%fitting_method == FittingMethodNO) &
      return 

    num_fit_atom = size(fitting%fitting_atom%idx)
    num_mov_atom = size(coord_fit(1,:))

    if (size(coord_fit(1,:)) /= num_mov_atom) &
        call error_msg('Run_Fitting> coord_ref, coord_fit : different size')
      
    coord_fit(1:3,1:num_mov_atom) = coord_mov(1:3,1:num_mov_atom)

    allocate(mass(size(coord_ref(1,:))))

    if (present(mass_ref)) then
      if (size(mass_ref) /= size(mass)) &
        call error_msg('Run_Fitting> coord_ref, mass_ref : different size')
      mass(:) = mass_ref(:)
    else
      mass(:) = 1.0_wp
    end if


    if (fitting%fitting_method == FittingMethodTR_ROT) then

      !  translation and rotation
      !
      call fit_trrot(num_fit_atom,             &
                     fitting%fitting_atom%idx, &
                     coord_ref,                &
                     mass,                     &
                     coord_mov,                &
                     fitting%rot_matrix,       &
                     fitting%com_ref,          &
                     fitting%com_mov,          &
                     fitting%rmsd,             &
                     fitting%ret_code)

    else if ((fitting%fitting_method == FittingMethodTR) .or. &
             (fitting%fitting_method == FittingMethodXYTR)) then 

      !  translation
      !
      call fit_trans(num_fit_atom,             &
                     fitting%fitting_atom%idx, &
                     coord_ref,                &
                     coord_mov,                &
                     mass,                     &
                     fitting%z_translation,    &
                     fitting%rot_matrix,       &
                     fitting%com_ref,          &
                     fitting%com_mov,          &
                     fitting%rmsd,             &
                     fitting%ret_code)

    else if (fitting%fitting_method == FittingMethodTR_ZROT) then

      !  translation + z-rotation
      !
      call fit_trans(num_fit_atom,             &
                     fitting%fitting_atom%idx, &
                     coord_ref,                &
                     coord_mov,                &
                     mass,                     &
                     fitting%z_translation,    &
                     fitting%rot_matrix,       &
                     fitting%com_ref,          &
                     fitting%com_mov,          &
                     fitting%rmsd,             &
                     fitting%ret_code)

      call fit_zrot(num_fit_atom,              &
                    fitting%fitting_atom%idx,  &
                    coord_ref,                 &
                    fitting%num_grids,         &
                    fitting%grid_size,         &
                    fitting%com_ref,           &
                    fitting%com_mov,           &
                    coord_mov,                 &
                    fitting%rot_matrix,        &
                    fitting%rmsd,              &
                    fitting%ret_code)


    else if (fitting%fitting_method == FittingMethodXYTR_ZROT) then

      !  translation + z-rotation
      !
      call fit_xytr_zrot2(num_fit_atom,        &
                     fitting%fitting_atom%idx, &
                     coord_ref,                &
                     coord_mov,                &
                     mass,                     &
                     fitting%rot_matrix,       &
                     fitting%com_ref,          &
                     fitting%com_mov,          &
                     fitting%rmsd,             &
                     fitting%ret_code)

    end if
 
    ! transform by translation & rotation information
    !
    call transform(num_mov_atom,            &
                   fitting%rot_matrix,      &
                   fitting%com_ref,         &
                   fitting%com_mov,         &
                   coord_fit)

    deallocate(mass)

    ! write(MsgOut,*) '                      RMSD = ', fitting%rmsd 
    ! write(MsgOut,*) ' '

    return

  end subroutine run_fitting
                         
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_rmsd
  !> @brief        output RMS deviation
  !! @authors      YS, NT
  !! @param[in]    file_unit_no : output destination
  !! @param[in]    structure_no : structure number
  !! @param[in]    fitting      : fitting information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_rmsd(file_unit_no, &
                      structure_no, &
                      fitting)

    ! formal arguments
    integer,                 intent(in)    :: file_unit_no
    integer,                 intent(in)    :: structure_no
    type(s_fitting),         intent(in)    :: fitting


    if (fitting%fitting_method == FittingMethodNO) &
      return

    write(file_unit_no, '(i10,1x,f8.3)') structure_no, fitting%rmsd

    return

  end subroutine out_rmsd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_trrot
  !> @brief        write a translational vector and a rotational matrix
  !! @authors      YS, NT, CK
  !! @param[in]    file_unit_no : unit number of out_trrot data
  !! @param[in]    structure_no : structure number
  !! @param[in]    fitting      : fitting information
  !! @note         
  !!    Translational vector
  !!      Tx  = cm2 - cm1
  !!
  !!    Rotational matrix (cf Goldstein)
  !!      Rxx =  cos(psi)*cos(phi)-cos(theta)*sin(phi)*sin(psi)
  !!      Rxy =  cos(psi)*sin(phi)+cos(theta)*cos(phi)*sin(psi)
  !!      Rxz =  sin(phi)*sin(theta)
  !!      Ryy = -sin(psi)*cos(phi)-cos(theta)*sin(phi)*cos(psi)
  !!      Ryx = -sin(psi)*cos(phi)+cos(theta)*cos(phi)*cos(psi)
  !!      Ryz =  cos(phi)*sin(theta)
  !!      Ryx =                    sin(theta)*sin(phi)
  !!      Ryy =                   -sin(theta)*cos(phi)
  !!      Rzz =                    cos(theta) 
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_trrot(file_unit_no, &
                       structure_no, &
                       fitting)

    ! formal arguments
    integer,                 intent(in) :: file_unit_no
    integer,                 intent(in) :: structure_no
    type(s_fitting),         intent(in) :: fitting

    ! local variables
    real(wp)                 :: rot_matrix(3,3), com_ref(3), com_mov(3), rotang


    if (fitting%fitting_method == FittingMethodNO) &
       return 

    rot_matrix(:,:) = fitting%rot_matrix(:,:)
    com_ref(:)      = fitting%com_ref(:)
    com_mov(:)      = fitting%com_mov(:)

    if ((fitting%fitting_method == FittingMethodTR_ZROT) .or. &
        (fitting%fitting_method == FittingMethodXYTR_ZROT)) then

       rotang = acos(rot_matrix(1,1))*RAD

       if (rot_matrix(2,1) < 0.0) rotang = -rotang

       write(file_unit_no, '(i10,1x,4(f8.3,1x))') &
                       structure_no,              &
                       com_mov(1) - com_ref(1),   &
                       com_mov(2) - com_ref(2),   &
                       com_mov(3) - com_ref(3),   &
                       rotang
    else
       write(file_unit_no, '(i10,1x,12(f8.3,1x))') &
                       structure_no,               &
                       com_mov(1) - com_ref(1),    &
                       com_mov(2) - com_ref(2),    &
                       com_mov(3) - com_ref(3),    &
                       rot_matrix(1,1), rot_matrix(1,2), &
                       rot_matrix(1,3), rot_matrix(2,1), &
                       rot_matrix(2,2), rot_matrix(2,3), &
                       rot_matrix(3,1), rot_matrix(3,2), &
                       rot_matrix(3,3)
    end if

    return

  end subroutine out_trrot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fit_trrot
  !> @brief        least square fitting of two conformations
  !! @authors      YS, NT, CK
  !! @param[in]    natom       : number of fitted atoms
  !! @param[in]    aidx        : fitted atom index
  !! @param[in]    refcoord    : reference coordinates
  !! @param[in]    mass_w      : mass weight of each atom
  !! @param[inout] movcoord    : moving coordinates
  !! @param[out]   rot_matrix  : rotation matrix for fitting
  !! @param[out]   com_ref     : center of mass of reference coordinate
  !! @param[out]   com_mov     : center of mass of moving coordinate
  !! @param[out]   rmsd        : root mean square deviation
  !! @param[out]   ierr        : return code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fit_trrot(natom, aidx, refcoord,                    &
                       mass_w, movcoord, rot_matrix,             &
                       com_ref, com_mov, rmsd,                   &
                       ierr)

    ! parameters
    integer,                 parameter     :: lwork = max(1,4*3-1)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: refcoord(:,:), mass_w(:)
    real(wp),                intent(in)    :: movcoord(:,:)
    integer,                 intent(in)    :: aidx(:)
    real(wp),                intent(out)   :: rot_matrix(3,3), com_ref(3)
    real(wp),                intent(out)   :: com_mov(3), rmsd
    integer,                 intent(out)   :: ierr

    ! local variables
    real(wp)                 :: total_mass, dref(1:3), dmov(1:3) 
    real(wp)                 :: dadd(1:3), dsub(1:3)
    real(wp)                 :: sym_matrix(4,4), eval(1:4), work(1:lwork)
    real(wp)                 :: evec(1:4)
    integer                  :: i, iatm


    if (natom == 0) then
      ierr = -1
      return 
    end if

    ! calculation mass of center
    !
    com_ref(1:3) = 0.0_wp
    com_mov(1:3) = 0.0_wp
    total_mass   = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      com_ref(1:3) = com_ref(1:3) + mass_w(iatm) * refcoord(1:3,iatm)
      com_mov(1:3) = com_mov(1:3) + mass_w(iatm) * movcoord(1:3,iatm)
      total_mass   = total_mass   + mass_w(iatm)
    end do

    if (total_mass <= 0.0_wp) then
      ierr = -2
      return
    end if

    total_mass   = 1.0_wp / total_mass
    com_ref(1:3) = com_ref(1:3) * total_mass
    com_mov(1:3) = com_mov(1:3) * total_mass

    ! calculation symmetric matrix
    !
    sym_matrix(1:4,1:4) = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      dref(1:3) = refcoord(1:3,iatm) - com_ref(1:3)
      dmov(1:3) = movcoord(1:3,iatm) - com_mov(1:3)
      dsub(1:3)   = mass_w(iatm) * (dref(1:3) - dmov(1:3))
      dadd(1:3)   = mass_w(iatm) * (dref(1:3) + dmov(1:3))
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

#ifdef LAPACK

    call dsyev('V', 'U', 4, sym_matrix, 4, eval, work, lwork, ierr)

#else

    call error_msg('Fit_Trrot> ERROR: This subroutine needs LAPACK.')

#endif

    if (ierr /= 0) then
      ierr = -2
      return
    end if

    ! calculation rotation matrix
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
    ! calculation of rmsd
    !
    rmsd = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      dref(1:3)  = refcoord(1:3,iatm) - com_ref(1:3)
      dmov(1:3)  = movcoord(1:3,iatm) - com_mov(1:3)
      dadd(1:3)    = rot_matrix(1:3,1)*dmov(1) + &
                     rot_matrix(1:3,2)*dmov(2) + &
                     rot_matrix(1:3,3)*dmov(3)
      rmsd = rmsd + (dref(1) - dadd(1))*(dref(1) - dadd(1)) + &
                    (dref(2) - dadd(2))*(dref(2) - dadd(2)) + &
                    (dref(3) - dadd(3))*(dref(3) - dadd(3))
    end do

    rmsd = sqrt(rmsd / real(natom,wp))

    ierr = 0
    return

  end subroutine fit_trrot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fit_trans
  !> @brief        fitting only by translation
  !! @authors      NT, CK
  !! @param[in]    natom      : number of fitted atoms
  !! @param[in]    aidx       : fitted atom index
  !! @param[in]    refcoord   : reference coordinates
  !! @param[inout] movcoord   : moving coordinates
  !! @param[in]    mass_w     : mass weight of each atom
  !! @param[in]    z_trans    : flag for z-translation of fitted coordinate
  !! @param[out]   rot_matrix : rotation matrix for fitting
  !! @param[out]   com_ref    : center of mass of reference coordinate
  !! @param[out]   com_mov    : center of mass of moving coordinate
  !! @param[out]   rmsd       : root mean square deviation
  !! @param[out]   ierr       : return code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fit_trans(natom, aidx, refcoord, &
                       movcoord, mass_w, z_trans, &
                       rot_matrix, com_ref, com_mov, &
                       rmsd,   ierr)

    ! formal arguments
    integer,                 intent(in)    :: natom
    integer,                 intent(in)    :: aidx(:)
    real(wp),                intent(in)    :: refcoord(:,:)
    real(wp),                intent(in)    :: movcoord(:,:)
    real(wp),                intent(in)    :: mass_w(:)
    logical,                 intent(in)    :: z_trans
    real(wp),                intent(out)   :: rot_matrix(3,3)
    real(wp),                intent(out)   :: com_ref(3)
    real(wp),                intent(out)   :: com_mov(3)
    real(wp),                intent(out)   :: rmsd
    integer,                 intent(out)   :: ierr

    ! local variables
    real(wp)                 :: total_mass, dref(1:3), dmov(1:3)
    integer                  :: i, iatm


    if (natom == 0) then
      ierr = -1
      return
    end if

    ! calculation mass of center
    !
    com_ref(1:3) = 0.0_wp
    com_mov(1:3) = 0.0_wp
    total_mass   = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      com_ref(1:3) = com_ref(1:3) + mass_w(iatm) * refcoord(1:3,iatm)
      com_mov(1:3) = com_mov(1:3) + mass_w(iatm) * movcoord(1:3,iatm)
      total_mass   = total_mass   + mass_w(iatm)
    end do

    if (total_mass <= 0.0_wp) then
      ierr = -2
      return
    end if

    total_mass   = 1.0_wp / total_mass
    com_ref(1:3) = com_ref(1:3) * total_mass
    com_mov(1:3) = com_mov(1:3) * total_mass

    ! calculation of rmsd
    !
    rot_matrix(1:3,1:3) = 0.0_wp
    rot_matrix(1,1)     = 1.0_wp
    rot_matrix(2,2)     = 1.0_wp
    rot_matrix(3,3)     = 1.0_wp

    rmsd = 0.0_wp

    if (.not. z_trans) then
      com_ref(3) = com_mov(3)
    end if

    do i = 1 , natom
      iatm = aidx(i)
      dref(1:3) = refcoord(1:3,iatm) - com_ref(1:3)
      dmov(1:3) = movcoord(1:3,iatm) - com_mov(1:3)
      rmsd = rmsd + (dref(1) - dmov(1))*(dref(1) - dmov(1)) + &
                    (dref(2) - dmov(2))*(dref(2) - dmov(2)) + &
                    (dref(3) - dmov(3))*(dref(3) - dmov(3))
    end do

    rmsd = sqrt(rmsd / real(natom,wp))

    ierr = 0
    return

  end subroutine fit_trans

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fit_zrot
  !> @brief        fitting only by z-rotation
  !! @authors      YS, NT, CK
  !! @param[in]    natom      : number of fitted atoms
  !! @param[in]    aidx       : fitted atom index
  !! @param[in]    refcoord   : reference coordinates
  !! @param[in]    num_grid   : number of grids
  !! @param[in]    grid_size  : grid size
  !! @param[in]    com_ref    : center of mass of reference coordinate
  !! @param[in]    com_mov    : center of mass of moving coordinate
  !! @param[inout] movcoord   : moving coordinates
  !! @param[out]   rot_matrix : rotation matrix for fitting
  !! @param[out]   rmsd       : root mean square deviation
  !! @param[out]   ierr       : return code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fit_zrot(natom, aidx, refcoord,                  &
                      num_grid, grid_size,                    &
                      com_ref, com_mov, movcoord, rot_matrix, &
                      rmsd,   ierr)
    
    ! formal arguments
    integer,                 intent(in)    :: natom
    integer,                 intent(in)    :: aidx(:)
    real(wp),                intent(in)    :: refcoord(:,:)
    integer,                 intent(in)    :: num_grid
    real(wp),                intent(in)    :: grid_size
    real(wp),                intent(in)    :: com_ref(3)
    real(wp),                intent(in)    :: com_mov(3)
    real(wp),                intent(in)    :: movcoord(:,:)
    real(wp),                intent(out)   :: rot_matrix(3,3)
    real(wp),                intent(out)   :: rmsd
    integer,                 intent(out)   :: ierr

    ! local variables
    real(wp)                 :: dref(1:3), dmov(1:3), rot_mov(1:3)
    real(wp)                 :: rotang, rmsdmin
    integer                  :: i, iatm, irot, irotmin


    ierr = 0

    rmsdmin = 10000000.0_wp

    do irot = -num_grid, num_grid

      rotang = grid_size * RAD * real(irot,wp)

      rot_matrix(1:3,1:3) = 0.0_wp
      rot_matrix(1,1)     = cos(rotang)
      rot_matrix(2,2)     = cos(rotang)
      rot_matrix(3,3)     = 1.0_wp
      rot_matrix(1,2)     =-sin(rotang)
      rot_matrix(2,1)     = sin(rotang)

      rmsd = 0.0_wp
      do i = 1, natom
        iatm = aidx(i)

        dref(1:3) = refcoord(1:3, iatm) - com_ref(1:3)
        dmov(1:3) = movcoord(1:3, iatm) - com_mov(1:3)

        rot_mov(1:3) = rot_matrix(1:3,1) * dmov(1) + &
                       rot_matrix(1:3,2) * dmov(2) + &
                       rot_matrix(1:3,3) * dmov(3)

        rmsd = rmsd + (dref(1) - rot_mov(1))*(dref(1) - rot_mov(1)) + &
                      (dref(2) - rot_mov(2))*(dref(2) - rot_mov(2)) + &
                      (dref(3) - rot_mov(3))*(dref(3) - rot_mov(3))
      end do
      rmsd = sqrt(rmsd / real(natom,wp))

      if (rmsd <= rmsdmin) then
        rmsdmin = rmsd
        irotmin = irot
      end if

    end do

    if ((irotmin  ==  num_grid) .or. (irotmin  ==  -num_grid)) &
      call error_msg('fitting> ERROR: more grids may be necessary!')

    rotang = grid_size * RAD * real(irotmin,wp)

    if (irotmin /= 0) then
      write(MsgOut,*) 'fit_zrot> angmin = ', grid_size*real(irotmin,wp)
      write(MsgOut,*) ' '
    end if

    rot_matrix(1:3,1:3) = 0.0_wp
    rot_matrix(1,1)     = cos(rotang)
    rot_matrix(2,2)     = cos(rotang)
    rot_matrix(3,3)     = 1.0_wp
    rot_matrix(1,2)     =-sin(rotang)
    rot_matrix(2,1)     = sin(rotang)

    rmsd = rmsdmin

    return

  end subroutine fit_zrot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    fit_xytr_zrot2
  !> @brief        least square fitting of two conformations
  !! @authors      CK
  !! @param[in]    natom       : number of fitted atoms
  !! @param[in]    aidx        : fitted atom index
  !! @param[in]    refcoord    : reference coordinates
  !! @param[in]    mass_w      : mass weight of each atom
  !! @param[in]    movcoord    : moving coordinates
  !! @param[out]   rot_matrix  : rotation matrix for fitting
  !! @param[out]   com_ref     : center of mass of reference coordinate
  !! @param[out]   com_mov     : center of mass of moving coordinate
  !! @param[out]   rmsd        : root mean square deviation
  !! @param[out]   ierr        : return code
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine fit_xytr_zrot2(natom, aidx, refcoord,               &
                       movcoord, mass_w, rot_matrix,             &
                       com_ref, com_mov, rmsd,                   &
                       ierr)

    ! parameters
    integer,                 parameter     :: lwork = max(1,5*2)

    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: refcoord(:,:)
    real(wp),                intent(in)    :: movcoord(:,:)
    real(wp),                intent(in)    :: mass_w(:)
    integer,                 intent(in)    :: aidx(:)
    real(wp),                intent(out)   :: rot_matrix(3,3), com_ref(3)
    real(wp),                intent(out)   :: com_mov(3), rmsd
    integer,                 intent(out)   :: ierr

    ! local variables
    real(wp)                 :: total_mass, dref(1:3), dmov(1:3) , dadd(1:3)
    real(wp)                 :: det1, det2
    real(wp)                 :: sym_matrix(2,2), u(2,2), vt(2,2)
    real(wp)                 :: eval(1:2), work(1:lwork)
    real(wp)                 :: evec(1:2)
    integer                  :: i, iatm


    if (natom == 0) then
      ierr = -1
      return 
    end if

    ! calculation mass of center
    !
    com_ref(1:3) = 0.0_wp
    com_mov(1:3) = 0.0_wp
    total_mass   = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      com_ref(1:3) = com_ref(1:3) + mass_w(iatm) * refcoord(1:3,iatm)
      com_mov(1:3) = com_mov(1:3) + mass_w(iatm) * movcoord(1:3,iatm)
      total_mass   = total_mass   + mass_w(iatm)
    end do

    if (total_mass <= 0.0_wp) then
      ierr = -2
      return
    end if

    total_mass   = 1.0_wp / total_mass
    com_ref(1:3) = com_ref(1:3) * total_mass
    com_mov(1:3) = com_mov(1:3) * total_mass

    ! calculation symmetric matrix
    !
    sym_matrix(1:2,1:2) = 0.0_wp
    rot_matrix(1:3,1:3) = 0.0_wp

    do i = 1, natom
      iatm = aidx(i)
      dref(1:2) = refcoord(1:2,iatm) - com_ref(1:2)
      dmov(1:2) = movcoord(1:2,iatm) - com_mov(1:2)

      sym_matrix(1, 1) = sym_matrix(1,1)+mass_w(iatm)*dref(1)*dmov(1)
      sym_matrix(1, 2) = sym_matrix(1,2)+mass_w(iatm)*dref(1)*dmov(2)
      sym_matrix(2, 1) = sym_matrix(2,1)+mass_w(iatm)*dref(2)*dmov(1)
      sym_matrix(2, 2) = sym_matrix(2,2)+mass_w(iatm)*dref(2)*dmov(2)

    end do

#ifdef LAPACK

    call dgesvd('A','A',2, 2, sym_matrix, 2, eval, U, 2, VT, 2, work, lwork, &
                ierr)
#else

    call error_msg('Fit_Trrot> ERROR: This subroutine needs LAPACK.')

#endif

    if (ierr /= 0) then
      ierr = -2
      return
    end if

    ! det(VT), det(U)
    det1 = (VT(1,1)*VT(2,2)-VT(1,2)*VT(2,1))
    det2 = (U(1,1)*U(2,2)-U(1,2)*U(2,1))
    if (det1*det2 < 0.0_wp) then
      eval(2) = -eval(2)
      U(1:2,2) = -U(1:2,2)
    end if
    rot_matrix(1,1) = U(1,1)*VT(1,1)+U(1,2)*VT(2,1)
    rot_matrix(2,1) = U(2,1)*VT(1,1)+U(2,2)*VT(2,1)
    rot_matrix(1,2) = U(1,1)*VT(1,2)+U(1,2)*VT(2,2)
    rot_matrix(2,2) = U(2,1)*VT(1,2)+U(2,2)*VT(2,2)
    rot_matrix(3,3) = 1.0_wp

    ! calculation of rmsd
    !
    rmsd = 0.0_wp

    com_ref(3) = com_mov(3)

    do i = 1, natom
      iatm = aidx(i)
      dref(1:3)  = refcoord(1:3,iatm) - com_ref(1:3)
      dmov(1:3)  = movcoord(1:3,iatm) - com_mov(1:3)
      dadd(1:3)    = rot_matrix(1:3,1)*dmov(1) + &
                     rot_matrix(1:3,2)*dmov(2) + &
                     rot_matrix(1:3,3)*dmov(3)
      rmsd = rmsd + (dref(1) - dadd(1))*(dref(1) - dadd(1)) + &
                    (dref(2) - dadd(2))*(dref(2) - dadd(2)) + &
                    (dref(3) - dadd(3))*(dref(3) - dadd(3))
    end do

    rmsd = sqrt(rmsd / real(natom,wp))

    ierr = 0
    return

  end subroutine fit_xytr_zrot2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    transform 
  !> @brief        rigid-body fitting
  !! @authors      YS, NT, CK
  !! @param[in]    natom      : number of atoms
  !! @param[in]    rot_matrix : rotation matrix
  !! @param[in]    com_ref    : center of mass of reference coordinate
  !! @param[in]    com_mov    : center of mass of moving coordinate
  !! @param[inout] coord      : coordinate data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
    
  subroutine transform(natom, rot_matrix, com_ref, com_mov, coord)
    
    ! formal arguments
    integer,                 intent(in)    :: natom
    real(wp),                intent(in)    :: rot_matrix(3,3)
    real(wp),                intent(in)    :: com_ref(3)
    real(wp),                intent(in)    :: com_mov(3)
    real(wp),                intent(inout) :: coord(:,:)

    ! local variables
    real(wp)                 :: dmov(1:3), rot_mov(1:3)
    integer                  :: iatm


    do iatm = 1, natom

      dmov(1:3)       = coord(1:3,iatm) - com_mov(1:3)
      rot_mov(1:3)    = rot_matrix(1:3,1) * dmov(1) + &
                        rot_matrix(1:3,2) * dmov(2) + &
                        rot_matrix(1:3,3) * dmov(3)
      coord(1:3,iatm) = rot_mov(1:3) + com_ref(1:3)

    end do

    return

  end subroutine transform

end module fitting_mod
