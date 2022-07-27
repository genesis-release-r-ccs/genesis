!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   ta_anlyze_mod
!> @brief   geometrical analysis of TM helices
!! @authors Daisuke Matsuoka (DM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module ta_analyze_mod

  use ta_option_mod
  use ta_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use fileio_mod
  use fileio_trj_mod
  use measure_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private 


  ! structure
  type, private :: s_helix
    type(s_selatoms)      :: atom
    real(wp)              :: axis(3)      ! unit vector parallel to the helix axis
  end type s_helix

  ! subroutines
  public  :: analyze
  private :: helix_axis
  private :: tilt_angle
  private :: dealloc_helix

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyses of a bilayer membrane and TM helices
  !! @authors      DM
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

 subroutine analyze(molecule,   &
                    trj_list,   &
                    trajectory, &
                    option,     &
                    output)

    ! formal arguments
    !
    type(s_molecule),   intent(inout) :: molecule
    type(s_trj_list),   intent(inout) :: trj_list
    type(s_trajectory), intent(inout) :: trajectory
    type(s_option),     intent(inout) :: option
    type(s_output),     intent(inout) :: output

    ! local variables
    !
    type(s_helix)       :: helix
    type(s_trj_file)    :: trj_in
    integer             :: nstru, irun, itrj, unit_no


    ! check-only
    if (option%check_only) &
      return

    call init_helix(option%TM_helix_atom, molecule, helix)

    ! open output files
    !
    call open_file(unit_no, output%outfile, IOFileOutputNew)

    nstru = 0

    do irun = 1, size(trj_list%md_steps)

      call open_trj(trj_in,                   &
                    trj_list%filenames(irun), &
                    trj_list%trj_format,      &
                    trj_list%trj_type,        &
                    IOFileInput)

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        call read_trj(trj_in, trajectory)

        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru          

          call helix_axis(trajectory%coord, helix)

           write(unit_no,'(I15, F15.3)') nstru, tilt_angle(helix)

        end if

      end do

      call close_trj(trj_in)

    end do

    call close_file(unit_no)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: Column 2: tilt angle of a helix to the z-axis (degree)'
    write(MsgOut,'(A)') ''


    call dealloc_helix(helix)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_helix
  !> @brief        extract only CA position of the TM helix
  !! @authors      DM
  !! @param[in]    selatoms : information about selected atom group
  !! @param[in]    molecule : molecule information
  !! @param[out]   helix    : helix information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_helix(selatoms, molecule, helix)

    ! formal arguments
    type(s_selatoms), intent(in)  :: selatoms
    type(s_molecule), intent(in)  :: molecule
    type(s_helix),    intent(out) :: helix

    ! local variable
    integer :: iatm, i, i_ca
    integer :: natm, n_ca

  
    if (allocated(helix%atom%idx)) &
      deallocate(helix%atom%idx)

    n_ca = 0
    natm = size(selatoms%idx)
    do i = 1, natm
      iatm = selatoms%idx(i)
      if (adjustl(molecule%atom_name(iatm)) == 'CA  ') &
        n_ca = n_ca + 1
    end do

    if (n_ca == 0) &
      call error_msg('helix_atom> No CA atoms are not found in the selection.')
  
    if (n_ca <= 3) &
      call error_msg('helix_atom> # of CA atoms are too small.')

    allocate(helix%atom%idx(n_ca))
  
    if (n_ca == natm) then
      helix%atom%idx = selatoms%idx

    else
      i_ca = 0
      do i = 1, natm
        iatm = selatoms%idx(i)
        if (adjustl(molecule%atom_name(iatm)) == 'CA  ') then
          i_ca = i_ca + 1
          helix%atom%idx(i_ca) = selatoms%idx(i)
        end if
      end do
    end if

    return

  end subroutine init_helix

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    helix_axis
  !> @brief        calculate a helix principal axis
  !! @authors      DM
  !! @param[in]    coord : atom coordinate
  !! @param[inout] helix : helix information
  !
  !! @note: Chotica et al. J. Mol. Biol. (1981) 145, 215
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine helix_axis(coord, helix)

    ! formal arguments
    real(wp),      intent(in)    :: coord(:,:)
    type(s_helix), intent(inout) :: helix

    ! local variables
    integer               :: nres, ires
    integer               :: i0, i1, i2
    integer               :: i, j

    real(wp), allocatable :: p(:,:)
    real(wp)              :: p_sum(3), p_ave(3)
    real(wp)              :: m(3,3)
    real(wp)              :: a(3)

    ! for LAPACK
    integer               :: INFO, LDA, LWORK, N
    character(LEN=1)      :: JOBZ, UPLO
    real(wp), allocatable :: WORK(:)
    real(wp)              :: eigvec(3)


    nres = size(helix%atom%idx)
    if (nres < 3) &
      call error_msg('helix_axis> # of the helix residues are less than 3.')

    ! define a circle in a plane normal to the helix axis
    ! p = r(i) + r(i+2) - 2 * r(i+1)
    !
    allocate(p(3, nres-2))

    do ires = 1, nres-2
      i0 = helix%atom%idx(ires)
      i1 = helix%atom%idx(ires+1)
      i2 = helix%atom%idx(ires+2)

      p(:,ires) = coord(:,i0) + coord(:,i2) - 2.0_wp * coord(:,i1)
    end do 

    !  N-2
    ! SIGMA[p(i) + p(i+2) - 2 * p(i+1)] = r(1) - r(2) - r(N-1) + r(N)
    !  i=1
    !
    p_sum(:) = (coord(:,helix%atom%idx(1))    - coord(:,helix%atom%idx(2))) + &
               (coord(:,helix%atom%idx(nres)) - coord(:,helix%atom%idx(nres-1)))
    p_ave(:) = p_sum(:) / real(nres-2, wp)

    ! inertia matrix
    !
    m = 0.0_wp
    do j = 1,3
      do i = 1, 3
        do ires = 1, nres-2
          m(i,j) = m(i,j) + (p(i,ires) - p_ave(i)) * (p(j,ires) - p_ave(j))
        end do
      end do
    end do

    ! calculate the vector parallel to the principal axis of the helix
    !
#ifdef LAPACK
    JOBZ  = 'V'
    UPLO  = 'L'
    N     = size(m,1)
    LDA   = size(m,1)
    LWORK = 3*N-1
    INFO  = 0

#ifdef _SINGLE

    allocate(WORK(LWORK))
    call ssyev(JOBZ, UPLO, N, m, LDA, eigvec, WORK, LWORK, INFO)

#else

    allocate(WORK(LWORK))
    call dsyev(JOBZ, UPLO, N, m, LDA, eigvec, WORK, LWORK, INFO)

#endif

    if (m(3,1) < 0.0_wp) then
      a(1:3) = -1.0_wp * m(1:3, 1)
    else
      a(1:3) = m(1:3, 1)
    end if

#else

    call error_msg('helix_axis> ERROR: This subroutine needs LAPACK')

#endif

    helix%axis(:) = a(:)

    ! deallocate memory
    !
    deallocate(WORK)
    deallocate(p)

    return

  end subroutine helix_axis

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    tilt_angle
  !> @brief        calculate a angle between a helix axis and the z-axis
  !! @authors      DM
  !! @param[in]    axis : a vector of a helix principal axis
  !
  !! @note: Chotica et al. J. Mol. Biol. (1981) 145, 215
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function tilt_angle(helix)

    ! formal arguments
    type(s_helix), intent(in) :: helix
    real(wp)                  :: tilt_angle

    tilt_angle = acos(helix%axis(3)) / RAD

    return

  end function tilt_angle

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    deallocate_helix
  !> @brief        deallocate helix information
  !! @authors      DM
  !! @param[inout] helix : helix information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine dealloc_helix(helix)

    ! format argument
    !
    type(s_helix), intent(inout) :: helix

    if (allocated(helix%atom%idx)) &
      deallocate(helix%atom%idx)

    return

    end subroutine dealloc_helix

end module ta_analyze_mod
