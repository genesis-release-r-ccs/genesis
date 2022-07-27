!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   lt_anlyze_mod
!> @brief   analysis of lipid bilayer thickness
!! @authors Daisuke Matsuoka (DM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module lt_analyze_mod

  use lt_option_mod
  use lt_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use molecules_str_mod
  use pbc_correct_mod
  use fileio_mod
  use fileio_trj_mod
  use measure_mod
  use messages_mod
  use constants_mod
 
  implicit none
  private 

  ! structure
  type, private :: s_membrane
    type(s_selatoms)      :: atom
    real(wp)              :: Zcrd_upper
    real(wp)              :: Zcrd_lower
  end type s_membrane

  ! subroutines
  public  :: analyze
  private :: bilayer
  private :: residue_com

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyses of a lipid bilayer thickness
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
    type(s_membrane)    :: membrane
    type(s_trj_file)    :: trj_in
    integer             :: nstru, irun, itrj
    integer             :: unit_no
    real(wp)            :: thickness


    ! check-only
    if (option%check_only) &
      return

    membrane%atom = option%membrane_atom

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

          call bilayer(molecule, trajectory, membrane)

          thickness = membrane%Zcrd_upper - membrane%Zcrd_lower
          write(unit_no, '(I15,3F10.3)') nstru, thickness,     &
                                         membrane%Zcrd_upper,  &
                                         membrane%Zcrd_lower
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
    write(MsgOut,'(A)') '    Column 2: lipid thickness (angstrom)'
    write(MsgOut,'(A)') '    Column 3: average z-coordinate of phosphorous atoms in the upper leaflet (angstrom)'
    write(MsgOut,'(A)') '    Column 4: average z-coordinate of phosphorous atoms in the lower leaflet (angstrom)'
    write(MsgOut,'(A)') ''


    return

  end subroutine analyze
 
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    bilayer
  !> @brief        move the center of the lipid bilayer to z=0
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] membrane   : membrane information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine bilayer(molecule, trajectory, membrane)

    type(s_molecule),        intent(in)    :: molecule
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_membrane),        intent(inout) :: membrane

    ! local variables
    integer                  :: i, j, k, iatm, jatm
    integer                  :: n_mematm, natm
    real(wp)                 :: com(3)
    real(wp)                 :: coord(3)
    real(wp)                 :: box_size_z

    integer                  :: nrs_pre, nmol
    integer                  :: first, last

    integer                  :: num_upper       ! number of lipid molecules in upper leaflet
    integer                  :: num_lower       ! number of lipid molecules in lower leaflet
    real(wp)                 :: mass_upper      ! total mass of selected group in upper leaflet
    real(wp)                 :: mass_lower      ! total mass of selected group in lower leaflet
    real(wp),    allocatable :: upper_crd(:,:)  ! coordinates of phosphorus atoms in upper leaflet
    real(wp),    allocatable :: lower_crd(:,:)  ! coordinates of phosphorus atoms in lower leaflet
    real(wp)                 :: upper_z         ! average value of z component of phosphorus atom
                                                ! in upper leaflet
    real(wp)                 :: lower_z         ! average value of z component of phosphorus atom
                                                ! in lower leaflet
    real(wp),    allocatable :: ul_idx(:)       ! phosphorus atom in upper leaflet => +1
                                                ! phosphorus atom in lower leaflet => -1
                                                ! other atoms                      =>  0
    integer                  :: idx_ul

    n_mematm = size(membrane%atom%idx)
    natm = molecule%num_atoms

    allocate(ul_idx(n_mematm))
    ul_idx(:) = 0.0_wp

    ! step 1: caluclate the center of mass of the membrane bilayer
    !         and move to the mid-plane of the membrane bilayer to z=0
    !
    com = compute_com(trajectory%coord, molecule%mass, membrane%atom%idx)

    do iatm = 1, natm
      coord(1:3) = trajectory%coord(1:3, iatm)
      trajectory%coord(1:3, iatm) = coord(1:3) - com(1:3)
    end do

    ! step 2: wrap molecules
    !
    call run_pbc_correct(PBCCModeMolecule, molecule, trajectory)


    ! step3: identify molecules in upper/lower leaflet
    !        use the COM position of each membrane moleucule for this purpose
    !
    nmol    = 0
    nrs_pre = -99999

    do i = 1, n_mematm
      iatm = membrane%atom%idx(i)

      if (molecule%residue_no(iatm) /= nrs_pre) then

        if (nmol > 0) then
          last  = i-1

          com = residue_com(trajectory%coord, molecule%mass, membrane%atom, first, last)

          ! calculation for the previous molecule
          !
          do j = first, last
            jatm = membrane%atom%idx(j)
            ul_idx(j) = sign(1.0_wp, com(3))
          end do
        end if   

        nmol = nmol + 1
        nrs_pre = molecule%residue_no(iatm)
        first = i

      else if (i == n_mematm) then
        last = i

        com = residue_com(trajectory%coord, molecule%mass, membrane%atom, first, last)

        ! calculation for the last molecule
        !
        do j = first, last
          jatm = membrane%atom%idx(j)
          ul_idx(j) = sign(1.0_wp, com(3))
        end do

      end if
    end do

    num_upper = count(ul_idx > 0.0_wp)
    num_lower = count(ul_idx < 0.0_wp)

    if (num_upper == 0) &
      call error_msg('BILAYER> There is no P atoms in upper leaflet')

    if (num_lower == 0) &
      call error_msg('BILAYER> There is no P atoms in lower leaflet')


#ifdef DEBUG
    write(MsgOut,*) nmol, num_upper, num_lower
#endif

    ! step4: measure bilayer thickness
    ! 
    allocate(upper_crd(3, num_upper))
    allocate(lower_crd(3, num_lower))

    j = 0
    k = 0
    mass_upper = 0.0_wp
    mass_lower = 0.0_wp

    do i = 1, n_mematm
      iatm = membrane%atom%idx(i)

      if (ul_idx(i) > 0.0_wp) then
        j = j + 1
        mass_upper = mass_upper + molecule%mass(iatm)
        upper_crd(1:3, j) = trajectory%coord(1:3, iatm) * molecule%mass(iatm)

      else if (ul_idx(i) < 0.0_wp) then
        k = k + 1
        mass_lower = mass_lower + molecule%mass(iatm)
        lower_crd(1:3, k) = trajectory%coord(1:3, iatm) * molecule%mass(iatm)

      end if
    end do

!     upper_z = sum(upper_crd(3,:)) / real(num_upper, wp)
!     lower_z = sum(lower_crd(3,:)) / real(num_lower, wp)
    upper_z = sum(upper_crd(3,:)) / mass_upper
    lower_z = sum(lower_crd(3,:)) / mass_lower

    membrane%Zcrd_upper = upper_z
    membrane%Zcrd_lower = lower_z

    ! step 5: deallocation
    !
    deallocate(ul_idx)
    deallocate(upper_crd, lower_crd)
  
    return

  end subroutine bilayer


  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    residue_com
  !> @brief        residue_com 
  !! @authors      DM
  !! @param[in]    coord    : atom coordinates
  !! @param[in]    mass     : atom masses
  !! @param[in]    selatoms : information about selected atoms
  !! @param[in]    first    : first atom position of the selected atom group
  !! @param[in]    last     : last atom position of the selected atom group
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function residue_com(coord, mass, selatoms, first, last)

    ! return value
    real(wp)                 :: residue_com(3)

    ! formal arguements
    real(wp),         intent(in)  :: coord(:,:)
    real(wp),         intent(in)  :: mass(:)
    type(s_selatoms), intent(in)  :: selatoms
    integer,          intent(in)  :: first, last

    ! local variables
    integer              :: i, natom
    integer, allocatable :: group(:) 


    if (allocated(group)) &
      deallocate(group)

    natom = last - first + 1
    allocate(group(natom))

    do i = 1, natom
      group(i) = selatoms%idx(first+i-1)
    end do
    residue_com = compute_com(coord, mass, group)

    deallocate(group)

    return

  end function residue_com

end module lt_analyze_mod
