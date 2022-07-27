!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   fa_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Daisuke Matsuoka (DM), Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module fa_analyze_mod

  use fa_option_str_mod
  use fileio_trj_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use fileio_mod
  use measure_mod
  use messages_mod
  use constants_mod
  use atom_libs_mod
 
  implicit none
  private

  ! structure
  type, private :: s_chromophore
    real(wp) :: com(3)      ! center of mass of chromophore
    real(wp) :: dipole(3)   ! transition dipole
  end type s_chromophore

  type, private :: s_fret
    real(wp) :: distance    ! dye-dye distance
    real(wp) :: kappa       ! orientation factor
    real(wp) :: efficiency
  end type s_fret

  ! subroutines
  public  :: analyze
  private :: analyze_fret
  private :: out_result

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories
  !! @authors      DM, NT
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    output     : output information
  !! @param[inout] option     : option information
  !! @param[inout] trajectory : trajectory information
  !! @note         PNAS (2005) 102, 2754-2759
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, output, option, trajectory)

    ! formal arguments
    type(s_molecule),   intent(in)    :: molecule
    type(s_trj_list),   intent(in)    :: trj_list
    type(s_output),     intent(in)    :: output
    type(s_option),     intent(inout) :: option
    type(s_trajectory), intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)          :: trj_in
    integer                   :: nstru, ifile, istep, num_trjfiles, fret_unit
    type(s_fret)              :: fret


    if (option%check_only) &
      return

    ! open output file
    !
    if (output%outfile /= '') &
      call open_file(fret_unit, output%outfile, IOFileOutputNew)

    ! analysis loop
    !
    nstru = 0
    num_trjfiles = size(trj_list%md_steps)

    do ifile = 1, num_trjfiles

      ! open trajectory file
      !
      call open_trj(trj_in, trj_list%filenames(ifile), &
                            trj_list%trj_format,       &
                            trj_list%trj_type, IOFileInput)

      do istep = 1, trj_list%md_steps(ifile)

        ! read trajectory
        !   coordinates of one MD snapshot are saved in trajectory%coord)
        !
        call read_trj(trj_in, trajectory)

        if (mod(istep, trj_list%ana_periods(ifile)) == 0) then

          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          call analyze_fret(trajectory, molecule, option, fret)
          call out_result (nstru, fret_unit, fret)

        end if

      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do


    ! close output file
    !
    call close_file(fret_unit)


    ! Output summary
    !
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
    write(MsgOut,'(A)') '    Column 1: Snapshot index'
    write(MsgOut,'(A)') '    Column 2: distance'
    write(MsgOut,'(A)') '    Column 3: orientation factor'
    write(MsgOut,'(A)') '    Column 4: FRET efficiency'
    write(MsgOut,'(A)') ''

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze_fret
  !> @brief        analysis distance RMS from near-native dimer configuration
  !! @authors      DM
  !! @param[in]    trajectory  : trajectory informatin
  !! @param[in]    molecule    : molecule informatin
  !! @param[in]    option      : option informatin
  !! @param[out]   fret        : FRET-related values
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze_fret(trajectory, molecule, option, fret)

    ! formal argments
    type(s_trajectory),   intent(in)    :: trajectory
    type(s_molecule),     intent(in)    :: molecule
    type(s_option),       intent(in)    :: option
    type(s_fret),         intent(out)   :: fret

    ! local variables
    type(s_chromophore)   :: chromophore(2)
    real(wp)              :: transfer_vec(3)
    real(wp)              :: dist
    real(wp)              :: diff(3)
    real(wp)              :: M_tensor(3,3)
    real(wp)              :: ratio
    real(wp)              :: denominator

    integer               :: i_chromo
    integer               :: i
    integer               :: iatm

    ! for LAPACK
    integer               :: INFO, LDA, LWORK, N
    character(len=1)      :: JOBZ, UPLO
    real(wp), allocatable :: WORK(:)
    real(wp)              :: eigvec(3)


    ! calculate the unit vector connecting the donor and the acceptor.
    !
    chromophore(1)%com(:) = compute_com(trajectory%coord,  &
                                        molecule%mass,     &
                                        option%analysis_atom(1)%idx)

    chromophore(2)%com(:) = compute_com(trajectory%coord,  &
                                        molecule%mass,     &
                                        option%analysis_atom(2)%idx)

    transfer_vec(:) = chromophore(1)%com(:) - chromophore(2)%com(:)
    dist = sqrt(dot_product(transfer_vec, transfer_vec))

    transfer_vec(:) = transfer_vec(:) / dist

    fret%distance = dist


    ! calculate transition dipole moment

    do i_chromo = 1, 2

      ! 1. compute the inertia moment
      !
      M_tensor(:,:) = 0.0_wp
    
      do i = 1, size(option%analysis_atom(i_chromo)%idx)
        iatm = option%analysis_atom(i_chromo)%idx(i)
        diff(:) = trajectory%coord(:, iatm) - chromophore(i_chromo)%com(:)

        M_tensor(1, 1) = M_tensor(1, 1) + (diff(2) * diff(2) + diff(3) * diff(3))
        M_tensor(1, 2) = M_tensor(1, 2) -  diff(1) * diff(2)
        M_tensor(1, 3) = M_tensor(1, 3) -  diff(1) * diff(3)

        M_tensor(2, 1) = M_tensor(2, 1) -  diff(2) * diff(1)
        M_tensor(2, 2) = M_tensor(2, 2) + (diff(3) * diff(3) + diff(1) * diff(1))
        M_tensor(2, 3) = M_tensor(2, 3) -  diff(2) * diff(3)

        M_tensor(3, 1) = M_tensor(3, 1) -  diff(3) * diff(1)
        M_tensor(3, 2) = M_tensor(3, 2) -  diff(3) * diff(2)
        M_tensor(3, 3) = M_tensor(3, 3) + (diff(1) * diff(1) + diff(2) * diff(2))
      end do

      ! 2. calculate principle axes
      !
#ifdef LAPACK
      JOBZ  = 'V'
      UPLO  = 'L'
      N     = size(M_tensor, 1)
      LDA   = size(M_tensor, 1)
      LWORK = 3 * N - 1
      INFO  = 0

      allocate(WORK(LWORK))

#ifdef _SINGLE
      call ssyev(JOBZ, UPLO, N, M_tensor, LDA, eigvec, WORK, LWORK, INFO)

#else
      call dsyev(JOBZ, UPLO, N, M_tensor, LDA, eigvec, WORK, LWORK, INFO)

#endif

#else
      call error_msg('FRET_analysis> ERROR: This subroutine requires LAPACK!')

#endif

      chromophore(i_chromo)%dipole(:) = M_tensor(:, 1)

      deallocate(WORK)
    end do

    ! compute oritentation factor
    !
    fret%kappa =   dot_product(chromophore(1)%dipole, chromophore(2)%dipole)  &
                 - 3.0_wp * dot_product(chromophore(1)%dipole, transfer_vec)  &
                          * dot_product(chromophore(2)%dipole, transfer_vec)

    ! FRET effieciency
    !
    ratio = fret%distance / option%Forster_radius

    denominator = (2.0_wp / 3.0_wp) * (1.0_wp / (fret%kappa ** 2)) * (ratio ** 6)

    fret%efficiency = 1.0_wp / (1.0_wp + denominator)

  end subroutine analyze_fret


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    out_result
  !> @brief        output analysis results
  !! @authors      NT
  !! @param[in]    sturct_no : sturcture number
  !! @param[in]    unit_no   : file unit number
  !! @param[in]    fret      : analysis results
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine out_result(struct_no, unit_no, fret)

    ! formal arguments
    integer,       intent(in)    :: struct_no
    integer,       intent(in)    :: unit_no
    type(s_fret),  intent(in) :: fret


    write(unit_no, '(i10,3(1x,f8.3))') struct_no, fret%distance,  &
                                                  fret%kappa,     &
                                                  fret%efficiency

    return

  end subroutine out_result

end module fa_analyze_mod
