!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fileio_gpr_mod
!> @brief   read and write GPR file
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module fileio_gpr_mod

  use fileio_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_gpr

    integer          :: type          = 0

    integer          :: num_pairs     = 0
    integer          :: num_bonds     = 0
    integer          :: num_angles    = 0
    integer          :: num_dihedrals = 0
    integer          :: num_impropers = 0

    ! bonds
    integer,          allocatable :: bond_list(:,:)
    real(wp),         allocatable :: bond_dist_min(:)
    real(wp),         allocatable :: bond_force_const(:)

    ! angles
    integer,          allocatable :: angl_list(:,:)
    real(wp),         allocatable :: angl_theta_min(:)
    real(wp),         allocatable :: angl_force_const(:)

    ! dihedrals
    integer,          allocatable :: dihe_list(:,:)
    real(wp),         allocatable :: dihe_phase(:)
    real(wp),         allocatable :: dihe_force_const(:)
    integer,          allocatable :: dihe_periodicity(:)

    ! impropers
    integer,          allocatable :: impr_list(:,:)
    real(wp),         allocatable :: impr_phase(:)
    real(wp),         allocatable :: impr_force_const(:)

    ! pairs
    integer,          allocatable :: pair_list(:,:)
    real(wp),         allocatable :: pair_eps(:)
    real(wp),         allocatable :: pair_dist_min(:)

    ! non-pairs
    real(wp)                      :: nonpair_eps      = 0.0_wp
    real(wp)                      :: nonpair_dist_min = 0.0_wp

  end type s_gpr

  ! parameters for allocatable variables
  integer,     public,  parameter :: GprBond = 1
  integer,     public,  parameter :: GprAngl = 2
  integer,     public,  parameter :: GprDihe = 3
  integer,     public,  parameter :: GprImpr = 4
  integer,     public,  parameter :: GprPair = 5

  ! subroutines
  public  :: input_gpr
  public  :: output_gpr
  public  :: init_gpr
  public  :: alloc_gpr
  public  :: dealloc_gpr
  public  :: dealloc_gpr_all

  private :: read_gpr
  private :: read_gpr_bond
  private :: read_gpr_angl
  private :: read_gpr_dihe
  private :: read_gpr_impr
  private :: read_gpr_pair
  private :: read_gpr_nonpair

  private :: write_gpr
  private :: write_gpr_bond
  private :: write_gpr_angl
  private :: write_gpr_dihe
  private :: write_gpr_impr
  private :: write_gpr_pair
  private :: write_gpr_nonpair


contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_gpr
  !> @brief        a driver subroutine for reading GPR file
  !! @authors      TM
  !! @param[in]    gpr_filename : filename of GPR file
  !! @param[out]   gpr          : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_gpr(gpr_filename, gpr)

    ! formal arguments
    character(*),            intent(in)    :: gpr_filename
    type(s_gpr),             intent(inout) :: gpr
    
    ! local variables
    integer                  :: file


    ! open GPR file
    !
    call open_file(file, gpr_filename, IOFileInput)

    ! read coordinate data from GPR file
    !
    call read_gpr(file, gpr)

    ! close GPR file
    !
    call close_file(file)

    return

  end subroutine input_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_gpr
  !> @brief        a driver subroutine for writing GPR file
  !! @authors      TM
  !! @param[in]    gpr_filename : filename of GPR file
  !! @param[in]    gpr          : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_gpr(gpr_filename, gpr)

    ! formal arguments
    character(*),            intent(in)    :: gpr_filename
    type(s_gpr),             intent(in)    :: gpr
    
    ! local variables
    integer                  :: file


    ! open GPR file
    !
    call open_file(file, gpr_filename, IOFileOutputNew)

    ! write coordinate data from GPR file
    !
    call write_gpr(file, gpr)

    ! close GPR file
    !
    call close_file(file)

    return

  end subroutine output_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_gpr
  !> @brief        initialize GPR information
  !! @authors      TM
  !! @param[out]   gpr : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_gpr(gpr)

    ! formal arguments
    type(s_gpr),             intent(inout) :: gpr


    gpr%num_bonds     = 0 
    gpr%num_angles    = 0
    gpr%num_dihedrals = 0
    gpr%num_impropers = 0
    gpr%num_pairs     = 0

    return

  end subroutine init_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    alloc_gpr
  !> @brief        allocate GPR information
  !! @authors      TM
  !! @param[inout] gpr      : structure of GPR information
  !! @param[in]    variable : selected variable
  !! @param[in]    var_size : size of the selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
 
  subroutine alloc_gpr(gpr, variable, var_size)

    ! formal arguments
    type(s_gpr),             intent(inout) :: gpr
    integer,                 intent(in)    :: variable
    integer,                 intent(in)    :: var_size

    ! local variables
    integer                  :: alloc_stat
    integer                  :: dealloc_stat


    ! initialize
    !
    alloc_stat   = 0
    dealloc_stat = 0

    ! allocate selected variables

    select case (variable)
      
    case (GprBond)

      if (allocated(gpr%bond_list)) then
        if (size(gpr%bond_list(1,:)) == var_size) &
          return
        deallocate(gpr%bond_list,        &
                   gpr%bond_dist_min,    & 
                   gpr%bond_force_const, &
                   stat = dealloc_stat)
      end if

      allocate(gpr%bond_list(2, var_size),     &
               gpr%bond_dist_min(var_size),    &
               gpr%bond_force_const(var_size), &
               stat = alloc_stat)

      gpr%bond_list  (1:2, 1:var_size) = 0
      gpr%bond_dist_min   (1:var_size) = 0.0_wp
      gpr%bond_force_const(1:var_size) = 0.0_wp

    case (GprAngl)

      if (allocated(gpr%angl_list)) then
        if (size(gpr%angl_list(1,:)) == var_size) &
          return
        deallocate(gpr%angl_list,        &
                   gpr%angl_theta_min,   &
                   gpr%angl_force_const, &
                   stat = dealloc_stat)
      end if

      allocate(gpr%angl_list(3, var_size),     &
               gpr%angl_theta_min(var_size),   &
               gpr%angl_force_const(var_size), &
               stat = alloc_stat)

      gpr%angl_list  (1:3, 1:var_size) = 0
      gpr%angl_theta_min  (1:var_size) = 0.0_wp
      gpr%angl_force_const(1:var_size) = 0.0_wp

    case (GprDihe)

      if (allocated(gpr%Dihe_list)) then
        if (size(gpr%dihe_list(1,:)) == var_size) &
          return
        deallocate(gpr%dihe_list,        &
                   gpr%dihe_phase,       &
                   gpr%dihe_force_const, &
                   gpr%dihe_periodicity, &
                   stat = dealloc_stat)
      end if

      allocate(gpr%dihe_list(4, var_size),     &
               gpr%dihe_phase(var_size),       &
               gpr%dihe_force_const(var_size), &
               gpr%dihe_periodicity(var_size), &
               stat = alloc_stat)

      gpr%dihe_list  (1:4, 1:var_size) = 0
      gpr%dihe_phase      (1:var_size) = 0.0_wp
      gpr%dihe_force_const(1:var_size) = 0.0_wp
      gpr%dihe_periodicity(1:var_size) = 0

    case (GprImpr)

      if (allocated(gpr%impr_list)) then
        if (size(gpr%impr_list(1,:)) == var_size) &
          return
        deallocate(gpr%impr_list,        &
                   gpr%impr_phase,       &
                   gpr%impr_force_const, &
                   stat = dealloc_stat)
      end if

      allocate(gpr%impr_list(4, var_size),     &
               gpr%impr_phase(var_size),       &
               gpr%impr_force_const(var_size), &
               stat = alloc_stat)

      gpr%impr_list  (1:4, 1:var_size) = 0
      gpr%impr_phase      (1:var_size) = 0.0_wp
      gpr%impr_force_const(1:var_size) = 0.0_wp

    case (GprPair)

      if (allocated(gpr%pair_list)) then
        if (size(gpr%pair_list(1,:)) == var_size) &
          return
        deallocate(gpr%pair_list,     &
                   gpr%pair_eps,      &
                   gpr%pair_dist_min, &
                   stat = dealloc_stat)
      end if

      allocate(gpr%pair_list(2, var_size),  &
               gpr%pair_eps(var_size),      &
               gpr%pair_dist_min(var_size), &
               stat = alloc_stat)

      gpr%pair_list(1:2, 1:var_size) = 0
      gpr%pair_eps      (1:var_size) = 0.0_wp
      gpr%pair_dist_min (1:var_size) = 0.0_wp

    case default

      call error_msg('Alloc_Gpr> bad variable')

    end select


    if (alloc_stat /= 0)   call error_msg_alloc
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine alloc_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_gpr
  !> @brief        deallocate GPR information
  !! @authors      TM
  !! @param[inout] gpr      : structure of GPR information
  !! @param[in]    variable : selected variable
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_gpr(gpr, variable)

    ! formal arguments
    type(s_gpr),             intent(inout) :: gpr
    integer,                 intent(in)    :: variable

    ! local variables
    integer                  :: dealloc_stat
   

    ! initialize
    !
    dealloc_stat = 0

    ! deallocate selected variable
    !
    select case (variable)

    case (GprBond)

      if (allocated(gpr%bond_list)) then
        deallocate( gpr%bond_list,        &
                    gpr%bond_dist_min,    &
                    gpr%bond_force_const, &
                    stat = dealloc_stat)
      end if

    case (GprAngl)

      if (allocated(gpr%angl_list)) then
        deallocate( gpr%angl_list,        &
                    gpr%angl_theta_min,   &
                    gpr%angl_force_const, &
                    stat = dealloc_stat)
      end if

    case (GprDihe)

      if (allocated(gpr%dihe_list)) then
        deallocate( gpr%dihe_list,        &
                    gpr%dihe_phase,       &
                    gpr%dihe_force_const, &
                    gpr%dihe_periodicity, &
                    stat = dealloc_stat)
      end if

    case (GprImpr)

      if (allocated(gpr%impr_list)) then
        deallocate( gpr%impr_list,        &
                    gpr%impr_phase,       &
                    gpr%impr_force_const, &
                    stat = dealloc_stat)
      end if

    case (GprPair)

      if (allocated(gpr%pair_list)) then
        deallocate( gpr%pair_list,        &
                    gpr%pair_eps,         &
                    gpr%pair_dist_min,    &
                    stat = dealloc_stat)
      end if

    case default

      call error_msg('Dealloc_Gpr> bad variable')

    end select

    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine dealloc_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    dealloc_gpr_all
  !> @brief        deallocate all GPR information
  !! @authors      TM
  !! @param[inout] gpr : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_gpr_all(gpr)

    ! formal arguments
    type(s_gpr),             intent(inout) :: gpr


    call dealloc_gpr(gpr, GprBond)
    call dealloc_gpr(gpr, GprAngl)
    call dealloc_gpr(gpr, GprDihe)
    call dealloc_gpr(gpr, GprImpr)
    call dealloc_gpr(gpr, GprPair)

    return

  end subroutine dealloc_gpr_all

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr
  !> @brief        read data from GPR file
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[out]   gpr  : gpr data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr(file, gpr)
       
    ! formal arugments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(inout) :: gpr
 
    ! local variables
    integer                  :: n1
    character(100)           :: key


    ! deallocate old data
    !
    call dealloc_gpr_all(gpr)
    call init_gpr(gpr)


    do while(.true.)

      read(file, '(I8,2X,A20)',err=100,end=100) n1, key

      if (key(1:4) == 'BOND') then
 
        ! read bond section
        !
        call read_gpr_bond(file, n1, gpr)
 
      else if (key(1:4) == 'ANGL') then
 
        ! read angle section
        !
        call read_gpr_angl(file, n1, gpr)
 
      else if (key(1:4) == 'DIHE') then
 
        ! read dihedral section
        !
         call read_gpr_dihe(file, n1, gpr)
 
      else if (key(1:4) =='IMPR') then
 
        ! read improper section
        !
        call read_gpr_impr(file, n1, gpr)

      else if (key(1:4) == 'PAIR') then

        ! read pair section
        !
        call read_gpr_pair(file, n1, gpr)

      else if (key(1:7) == 'NONPAIR') then

        ! read non-pair section
        !
        call read_gpr_nonpair(file, gpr)
 
      end if

    end do
100 continue


    ! write summary of PSF information
    !
    if (main_rank) then
      
      write(MsgOut,'(A)') 'Read_Gpr> Summary of GPR file'

      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_bonds       = ', gpr%num_bonds,     &
           '  num_angles      = ', gpr%num_angles
      write(MsgOut,'(A20,I10,A20,I10)')               &
           '  num_dihedrals   = ', gpr%num_dihedrals, &
           '  num_impropers   = ', gpr%num_impropers
      write(MsgOut,'(A20,I10)')                       &
           '  num_pairs       = ', gpr%num_pairs
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_bond
  !> @brief        read bonds
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    nbnd : number of bonds
  !! @param[inout] gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr_bond(file, nbnd, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nbnd
    type(s_gpr),             intent(inout) :: gpr

    ! local variables
    integer                  :: i


    ! Allocate Bond variables
    !
    call alloc_gpr(gpr, GprBond, nbnd)


    ! Read Bond
    !
    gpr%num_bonds = nbnd

    read(file, '(2i6,2x,2e13.5)') (gpr%bond_list(1, i),     &
                                   gpr%bond_list(2, i),     &
                                   gpr%bond_dist_min(i),    &
                                   gpr%bond_force_const(i), &
                                   i = 1, nbnd)

    return

  end subroutine read_gpr_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_angl
  !> @brief        read angles
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    nang : number of angles
  !! @param[inout] gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr_angl(file, nang, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nang
    type(s_gpr),             intent(inout) :: gpr

    ! local variables
    integer                  :: i


    ! Allocate Angle variables
    !
    call alloc_gpr(gpr, GprAngl, nang)


    ! Read Angl
    !
    gpr%num_angles = nang

    read(file, '(3i6,2x,2e13.5)') (gpr%angl_list(1, i),     &
                                   gpr%angl_list(2, i),     &
                                   gpr%angl_list(3, i),     &
                                   gpr%angl_theta_min(i),   &
                                   gpr%angl_force_const(i), &
                                   i = 1, nang)
    
    return

  end subroutine read_gpr_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_dihe
  !> @brief        read dihedrals
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    ndih : number of dihedrals
  !! @param[inout] gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr_dihe(file, ndih, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: ndih
    type(s_gpr),             intent(inout) :: gpr

    ! local variables
    integer                  :: i


    ! Allocate Dihedral variables
    !
    call alloc_gpr(gpr, GprDihe, ndih)


    ! read Dihedral
    !
    gpr%num_dihedrals = ndih

    read(file, '(4i6,2x,2e13.5,i2)') (gpr%dihe_list(1, i),     &
                                      gpr%dihe_list(2, i),     &
                                      gpr%dihe_list(3, i),     &
                                      gpr%dihe_list(4, i),     &
                                      gpr%dihe_phase(i),       &
                                      gpr%dihe_force_const(i), &
                                      gpr%dihe_periodicity(i), &
                                      i = 1, ndih)

    return

  end subroutine read_gpr_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_impr
  !> @brief        read impropers
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    nimp : number of impropers
  !! @param[inout] gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr_impr(file, nimp, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: nimp
    type(s_gpr),             intent(inout) :: gpr

    ! local variables
    integer                  :: i


    ! Allocate Improper Dihedral variables
    !
    call alloc_gpr(gpr, GprImpr, nimp)


    ! read Improper
    !
    gpr%num_impropers = nimp

    read(file, '(4i6,2x,2e13.5)') (gpr%impr_list(1, i),     &
                                   gpr%impr_list(2, i),     &
                                   gpr%impr_list(3, i),     &
                                   gpr%impr_list(4, i),     &
                                   gpr%impr_phase(i),       &
                                   gpr%impr_force_const(i), &
                                   i = 1, nimp)


    return

  end subroutine read_gpr_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_pair
  !> @brief        read pairlist
  !! @authors      TM
  !! @param[in]    file  : unit number of GPR file
  !! @param[in]    npair : number of pair lists
  !! @param[inout] gpr   : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_gpr_pair(file, npair, gpr)
  
    ! formal arguments
    integer,                 intent(in)    :: file
    integer,                 intent(in)    :: npair
    type(s_gpr),             intent(inout) :: gpr

    ! local variables
    integer                  :: i


    ! Allocate Pair variables
    !
    call alloc_gpr(gpr, GprPair, npair)


    ! Read Pairlist
    !
    gpr%num_pairs = npair

    read(file, '(2i6,2x,2e13.5)') (gpr%pair_list(1, i),  &
                                   gpr%pair_list(2, i),  &
                                   gpr%pair_eps(i),      &
                                   gpr%pair_dist_min(i), &
                                   i = 1, npair)

    return

  end subroutine read_gpr_pair

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_gpr_nonpair
  !> @brief        read non-pairlist
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[inout] gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_gpr_nonpair(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(inout) :: gpr


    read(file, '(2e13.5)') gpr%nonpair_eps, gpr%nonpair_dist_min

    return

  end subroutine read_gpr_nonpair

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr
  !> @brief        write data to GPR file
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : gpr data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr


    ! write summary of GPR information
    !
    write(MsgOut,*) 'Write_Gpr> Number of Bonds         = ', gpr%num_bonds
    write(MsgOut,*) 'Write_Gpr> Number of Angles        = ', gpr%num_angles
    write(MsgOut,*) 'Write_Gpr> Number of Dihedrals     = ', gpr%num_dihedrals
    write(MsgOut,*) 'Write_Gpr> Number of Impropers     = ', gpr%num_impropers
    write(MsgOut,*) 'Write_Gpr> Number of Pairs         = ', gpr%num_pairs
    write(MsgOut,*) ' '


    !  write bond section
    !
    call write_gpr_bond(file, gpr)

    !  write angle section
    !
    call write_gpr_angl(file, gpr)

    !  write dihedral section
    !
    call write_gpr_dihe(file, gpr)

    !  write improper section
    !
    call write_gpr_impr(file, gpr)

    !  write pair section
    !
    call write_gpr_pair(file, gpr)

    !  write nonpair section
    !
    call write_gpr_nonpair(file, gpr)


    return

  end subroutine write_gpr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_bond
  !> @brief        write bonds
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr_bond(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr

    ! local variables
    integer                  :: i, nbnd


    ! write the number of bonds
    !
    if (allocated(gpr%bond_list)) then
      nbnd = size(gpr%bond_list(1,:))
    else
      nbnd = 0
    end if

    write(file, '(i8,a)') nbnd, ' !BOND'
    

    ! write bond
    !
    write(file, '(2i6,2x,2e13.5)') (gpr%bond_list(1, i),     &
                                    gpr%bond_list(2, i),     &
                                    gpr%bond_dist_min(i),    &
                                    gpr%bond_force_const(i), &
                                    i = 1, nbnd)
    write(file,'(a)') ' '

    return

  end subroutine write_gpr_bond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_angl
  !> @brief        write angles
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr_angl(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr

    ! local variables
    integer                  :: i, nang


    ! write the number of angles
    !
    if (allocated(gpr%angl_list)) then
      nang = size(gpr%angl_list(1,:))
    else
      nang = 0
    end if

    write(file, '(i8,a)') nang, ' !ANGL'
    

    ! write angle
    !
    write(file, '(3i6,2x,2e13.5)') (gpr%angl_list(1, i),     &
                                    gpr%angl_list(2, i),     &
                                    gpr%angl_list(3, i),     &
                                    gpr%angl_theta_min(i),   &
                                    gpr%angl_force_const(i), &
                                    i = 1, nang)
    write(file,'(a)') ' '

    return

  end subroutine write_gpr_angl

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_dihe
  !> @brief        write dihedrals
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr_dihe(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr

    ! local variables
    integer                  :: i, ndih


    ! write the number of angles
    !
    if (allocated(gpr%dihe_list)) then
      ndih = size(gpr%dihe_list(1,:))
    else
      ndih = 0
    end if

    write(file, '(i8,a)') ndih, ' !DIHE'
    

    ! write dihe
    !
    write(file, '(4i6,2x,2e13.5,i2)') (gpr%dihe_list(1, i),     &
                                       gpr%dihe_list(2, i),     &
                                       gpr%dihe_list(3, i),     &
                                       gpr%dihe_list(4, i),     &
                                       gpr%dihe_phase(i),       &
                                       gpr%dihe_force_const(i), &
                                       gpr%dihe_periodicity(i), &
                                       i = 1, ndih)
    write(file,'(a)') ' '


    return

  end subroutine write_gpr_dihe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_impr
  !> @brief        write impropers
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr_impr(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr

    ! local variables
    integer                  :: i, nimp


    ! write the number of impropers
    !
    if (allocated(gpr%impr_list)) then
      nimp = size(gpr%impr_list(1,:))
    else
      nimp = 0
    end if

    write(file, '(i8,a)') nimp, ' !IMPR'
    

    ! Write impr
    !
    write(file, '(4i6,2x,2e13.5)') (gpr%impr_list(1, i),     &
                                    gpr%impr_list(2, i),     &
                                    gpr%impr_list(3, i),     &
                                    gpr%impr_list(4, i),     &
                                    gpr%impr_phase(i),       &
                                    gpr%impr_force_const(i), &
                                    i = 1, nimp)

    write(file,'(a)') ' '

    return

  end subroutine write_gpr_impr

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_pair
  !> @brief        write pairlist
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine write_gpr_pair(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr

    ! local variables
    integer                  :: i, npair


    ! write the number of pairs
    !
    if (allocated(gpr%pair_list)) then
      npair = size(gpr%pair_list(1,:))
    else
      npair = 0
    end if

    write(file, '(i8,a)') npair, ' !PAIR'


    ! write pairs
    !
    write(file, '(2i6,2x,2e13.5)') (gpr%pair_list(1, i),  &
                                    gpr%pair_list(2, i),  &
                                    gpr%pair_eps(i),      &
                                    gpr%pair_dist_min(i), &
                                    i = 1, npair)
    write(file,'(a)') ' '

    return

  end subroutine write_gpr_pair

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    write_gpr_nonpair
  !> @brief        write non-pairlist
  !! @authors      TM
  !! @param[in]    file : unit number of GPR file
  !! @param[in]    gpr  : structure of GPR information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine write_gpr_nonpair(file, gpr)

    ! formal arguments
    integer,                 intent(in)    :: file
    type(s_gpr),             intent(in)    :: gpr


    write(file, '(8x,a)') ' !NONPAIR'
    write(file, '(2e13.5)') gpr%nonpair_eps, gpr%nonpair_dist_min
    write(file,'(a)') ' '

    return
  
  end subroutine write_gpr_nonpair

end module fileio_gpr_mod
