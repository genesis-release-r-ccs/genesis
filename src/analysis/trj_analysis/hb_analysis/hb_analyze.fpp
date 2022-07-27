!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hb_analyze_mod
!> @brief   run analyzing trajectories
!! @authors Daisuke Matsuoka (DM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hb_analyze_mod

  use hb_option_str_mod
  use fileio_trj_mod
  use measure_mod
  use trajectory_str_mod
  use output_str_mod
  use molecules_str_mod
  use select_atoms_mod
  use select_atoms_str_mod
  use fileio_mod
  use messages_mod
  use constants_mod
  use atom_libs_mod
 
  implicit none
  private

  ! structure
  type, private :: s_polar_atom
    integer :: atom_no
    integer :: num_hydrogen
    integer :: hydrogen_atom(4)   ! e.g. NH4+ has four 'H-bond' donor.
  end type s_polar_atom

  type, private :: s_partner
    logical :: solvent
    integer :: atom_no
  end type s_partner

  type, private :: s_hb_info
    integer  :: partner
    real(wp) :: hb_dist
    real(wp) :: dha_angle
    real(wp) :: hda_angle
  end type s_hb_info

  ! subroutines
  public  :: analyze
  private :: get_polar_atom
  private :: examine_hbond
  private :: setup_hb_partner_list
  private :: is_polar_atom
  private :: print_output_info

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        run analyzing trajectories (analyze direct H-bonds)
  !! @authors      DM
  !! @param[in]    molecule   : molecule information
  !! @param[in]    trj_list   : trajectory file list information
  !! @param[in]    option     : option information
  !! @param[in]    output     : output information
  !! @param[inout] trajectory : trajectory information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine analyze(molecule, trj_list, trajectory, option, output)

    ! formal arguments
    type(s_molecule), target, intent(in)    :: molecule
    type(s_trj_list),         intent(in)    :: trj_list
    type(s_option),           intent(in)    :: option
    type(s_output),           intent(in)    :: output
    type(s_trajectory),       intent(inout) :: trajectory

    ! local variables
    type(s_trj_file)          :: trj_in
    integer                   :: nstru, ifile, istep, num_trjfiles
    integer                   :: iatm, jatm, idx, jdx
    integer                   :: out_unit, hb_out

    integer,          pointer :: numa(:)
    integer,          pointer :: numr(:)
    character(len=4), pointer :: nama(:)
    character(len=6), pointer :: namr(:)
    character(len=4), pointer :: seg(:)

    type(s_polar_atom), allocatable :: target_group(:)
    type(s_polar_atom), allocatable :: analysis_group(:)
    integer,            allocatable :: partner_list(:)
    type(s_partner),    allocatable :: partner_atom(:)
    integer,            allocatable :: hb_count(:,:)
    integer                         :: partner_idx
    logical                         :: hbond
    integer                         :: hb_total
    integer,            allocatable :: continue_hbond(:,:)
    type(s_hb_info)                 :: hb_list
    real(wp)                        :: dist, dha_angle, hda_angle


    ! formats
100 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',a4,1x,a6,1x,i7,1x,a4' | ',F6.3,2F9.3)
102 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',a4,1x,a6,1x,i7,1x,a4' | ',i10)
104 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',a4,1x,a6)
106 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',a4,1x,a6,1x,i7,1x,a4)
108 format('snapshot',i10,' : ',i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',a4,1x,a6,1x,i7,1x,a4)


    numa => molecule%atom_no
    nama => molecule%atom_name
    namr => molecule%residue_name
    numr => molecule%residue_no
    seg  => molecule%segment_name


    ! check option
    !
    if (option%boundary_type == BoundaryTypePBC) then
      if (trj_list%trj_type /= TrjTypeCoorBox) then
        call error_msg('hb_analysis> when boundary_type = PBC, TrjType must be COOR+BOX.')
      end if
    end if

    if (output%hb_listfile /= '') then
      call open_file(hb_out, output%hb_listfile, IOFileOutputNew)
    end if
    call open_file(out_unit, output%outfile, IOFileOutputNew)


    ! setup polar atoms (O, N)
    !

    ! for target_atom
    call get_polar_atom(molecule, option%target_atom, target_group)
    call setup_hb_partner_list(molecule, option, target_group,  &
                               partner_list, partner_atom)

#ifdef DEBUG
    write(MsgOut,'(A)') 'HB_Analyze> polar atoms are extracted from the target_atom group.'

    do iatm = 1, size(target_group)
      write(MsgOut,'(I7,$)') target_group(iatm)%atom_no
      if (mod(iatm, 10) == 0 .and. iatm /= size(target_group)) then
        write(MsgOut,'(A)') ''
      end if
    end do
    write(MsgOut,'(A/)') ''
#endif

    ! for analysis_atom
    call get_polar_atom(molecule, option%analysis_atom, analysis_group)

#ifdef DEBUG
    write(MsgOut,'(A)') 'HB_Analyze> polar atoms are extracted from the analysis_atom group.'

    do iatm = 1, size(analysis_group)
      write(MsgOut,'(I7,$)') analysis_group(iatm)%atom_no
      if (mod(iatm, 10) == 0 .and. iatm /= size(analysis_group)) then
        write(MsgOut,'(A)') ''
      end if
    end do
    write(MsgOut,'(A/)') ''
#endif

    if (option%check_only) &
      return

    ! analysis loop
    !
    nstru = 0
    num_trjfiles = size(trj_list%md_steps)

    allocate(hb_count(size(partner_atom), size(analysis_group)))
    allocate(continue_hbond(size(target_group), size(analysis_group)))

    hb_count(:,:) = 0
    continue_hbond(:,:) = 0

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

          hb_total = 0
          nstru = nstru + 1
          write(MsgOut,*) '      number of structures = ', nstru

          do idx = 1, size(analysis_group)
            do jdx = 1, size(target_group)

              call examine_hbond(analysis_group(idx), &
                                 target_group(jdx),   &
                                 trajectory,          &
                                 option,              &
                                 hbond, hb_list)

              if (output%hb_listfile /= '') then
                if (hbond) then
                  iatm = analysis_group(idx)%atom_no
                  jatm = target_group(jdx)%atom_no

                  write(hb_out, 100) &
                       nstru, &
                       nama(iatm), namr(iatm), numr(iatm), seg(iatm), &
                       nama(jatm), namr(jatm), numr(jatm), seg(jatm), &
                       hb_list%hb_dist, hb_list%dha_angle, hb_list%hda_angle
                end if
              end if

              select case (option%output_type)

                case (HBOutputModeCountSnap)
                  if (hbond) then
                    hb_total = hb_total + 1
                  end if

                case (HBOutputModeCountAtom)
                  if (hbond) then
                    partner_idx = partner_list(jdx)
                    hb_count(partner_idx, idx) = hb_count(partner_idx, idx) + 1
                  end if

                case (HBOutputModeLifetime)
                  if (hbond) then
                    continue_Hbond(jdx, idx) = continue_Hbond(jdx, idx) + 1

                  else
                    if (continue_Hbond(jdx, idx) > 0) then

                      iatm = analysis_group(idx)%atom_no
                      jatm = target_group(jdx)%atom_no

                      write(out_unit,108) &
                           nstru - 1, continue_Hbond(jdx, idx), &
                           nama(iatm), namr(iatm), numr(iatm), seg(iatm), &
                           nama(jatm), namr(jatm), numr(jatm), seg(jatm)

                      continue_Hbond(jdx, idx) = 0
                    end if
                  end if
              end select  

            end do
          end do

          if (option%output_type == HBOutputModeCountSnap) then
            write(out_unit,'(i10,2x,i10)') nstru, hb_total
          end if

        end if
      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    ! output the results
    !
    select case (option%output_type)

    case (HBOutputModeCountAtom)
      do idx = 1, size(analysis_group)
        do partner_idx = 1, size(partner_atom)

          if (hb_count(partner_idx, idx) == 0)  &
            cycle

          iatm = analysis_group(idx)%atom_no
          jatm = partner_atom(partner_idx)%atom_no

          if (partner_atom(partner_idx)%solvent) then
            write(out_unit, 104) hb_count(partner_idx, idx), &
                                 nama(iatm), namr(iatm), numr(iatm), seg(iatm), &
                                 nama(jatm), namr(jatm)

          else
            write(out_unit, 106) hb_count(partner_idx, idx), &
                                 nama(iatm), namr(iatm), numr(iatm), seg(iatm), &
                                 nama(jatm), namr(jatm), numr(jatm), seg(jatm)
          end if
        end do
      end do

    case (HBOutputModeLifetime)
      do idx = 1, size(analysis_group)
        do jdx = 1, size(target_group)

          iatm = analysis_group(idx)%atom_no
          jatm =   target_group(jdx)%atom_no

          if (continue_Hbond(jdx, idx) > 0) &
            write(out_unit, 108) nstru, continue_Hbond(jdx, idx), &
                                 nama(iatm), namr(iatm), numr(iatm), seg(iatm), &
                                 nama(jatm), namr(jatm), numr(jatm), seg(jatm)
        end do
      end do
    end select

    ! close output file
    !
    call close_file(out_unit)

    ! Output summary
    !
    call print_output_info(output, option)

    ! deacllocate
    !
    nullify(numa, numr, namr, nama, seg)
    deallocate(analysis_group, target_group)
    deallocate(hb_count, continue_hbond)
    deallocate(partner_list, partner_atom)

    return

  end subroutine analyze

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_polar_atoms
  !> @brief        extract only polar atoms from the specified atom group
  !! @authors      DM
  !! @param[in]    molecule    : molecule information
  !! @param[in]    select_atom : selected atom group
  !! @param[out]   polar_atom  : hydrogen bond-able atom information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
              
  subroutine get_polar_atom(molecule, select_atom, polar_atom)

    ! formal variables
    type(s_molecule),                intent(in)  :: molecule
    type(s_selatoms),                intent(in)  :: select_atom
    type(s_polar_atom), allocatable, intent(out) :: polar_atom(:)

    ! local variables
    integer              :: num_polar
    integer              :: iatm, natom
    integer              :: idx, atom_no
    integer, allocatable :: atom_to_idx(:)

    real(wp)             :: mass, mass1, mass2
    integer              :: ibond, num_hydrogen
    integer              :: bond_atom1, bond_atom2
    character(6)         :: atomcls, atomcls1, atomcls2


    ! extract polar atoms only from the selected atom group
    !
    num_polar = 0
    natom = size(select_atom%idx)

    do iatm = 1, natom
      atom_no = select_atom%idx(iatm)
      atomcls = molecule%atom_cls_name(atom_no)
      mass    = molecule%mass(atom_no)

      if (is_polar_atom(mass, atomcls)) then
        num_polar = num_polar + 1
      end if
    end do

    if (allocated(polar_atom)) &
      deallocate(polar_atom)

    allocate(polar_atom(num_polar))
    polar_atom(:)%num_hydrogen = 0

    allocate(atom_to_idx(molecule%num_atoms))
    atom_to_idx(:) = -99999

    num_polar = 0
    do iatm = 1, natom
      atom_no = select_atom%idx(iatm)
      atomcls = molecule%atom_cls_name(atom_no)
      mass    = molecule%mass(atom_no)

      if (is_polar_atom(mass, atomcls)) then

        num_polar = num_polar + 1
        polar_atom(num_polar)%atom_no = atom_no
        atom_to_idx(atom_no) = num_polar

      end if
    end do

    ! get hydrogen atoms bonded to the polar atoms
    !
    do ibond = 1, molecule%num_bonds
      bond_atom1 = molecule%bond_list(1, ibond)
      atomcls1   = molecule%atom_cls_name(bond_atom1)
      mass1      = molecule%mass(bond_atom1)

      bond_atom2 = molecule%bond_list(2, ibond)
      atomcls2   = molecule%atom_cls_name(bond_atom2)
      mass2      = molecule%mass(bond_atom2)

      if (atomic_number(mass1, atomcls1(1:6)) /= 1 .and. &
          atomic_number(mass2, atomcls2(1:6)) /= 1)      &
        cycle

      if (atom_to_idx(bond_atom1) > 0) then
        idx = atom_to_idx(bond_atom1)
        polar_atom(idx)%num_hydrogen = polar_atom(idx)%num_hydrogen + 1

        iatm = polar_atom(idx)%num_hydrogen
        polar_atom(idx)%hydrogen_atom(iatm) = bond_atom2

      else if (atom_to_idx(bond_atom2) > 0) then
        idx = atom_to_idx(bond_atom2)
        polar_atom(idx)%num_hydrogen = polar_atom(idx)%num_hydrogen + 1

        iatm = polar_atom(idx)%num_hydrogen
        polar_atom(idx)%hydrogen_atom(iatm) = bond_atom1
      end if

    end do

#ifdef DEBUB
    do idx = 1, size(polar_atom)
      atom_no      = polar_atom(idx)%atom_no
      num_hydrogen = polar_atom(idx)%num_hydrogen
      write(MsgOut,'(i7,i2,5A6)') idx, num_hydrogen,             &
                                  molecule%atom_name(atom_no),   &
                                 (molecule%atom_name(polar_atom(idx)%hydrogen_atom(iatm)), &
                                                     iatm=1, num_hydrogen)
    end do
#endif

    deallocate(atom_to_idx)

    return

  end subroutine get_polar_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hb_partner_list
  !> @brief        make H-bond partner list
  !! @authors      DM
  !! @param[in]    molecule     : molecule information
  !! @param[in]    option       : option information
  !! @param[in]    polar_atom   : H-bondable atom information
  !! @param[out]   parnter_idx  : polar atom index to partner index
  !! @param[out]   partner_atom : H-bond partner information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hb_partner_list(molecule, option, polar_atom, &
                                   partner_idx, partner_atom)

    ! formal arguments
    type(s_molecule),                intent(in)    :: molecule
    type(s_option),                  intent(in)    :: option
    type(s_polar_atom),              intent(in)    :: polar_atom(:)
    integer,            allocatable, intent(out)   :: partner_idx(:)
    type(s_partner),    allocatable, intent(out)   :: partner_atom(:)

    ! local variable
    integer                      :: n_solvent, i_solvent
    integer,         allocatable :: n_polar(:)
    integer                      :: i_polar
    logical,         allocatable :: detect_solvent(:)
    integer,         allocatable :: solvent_polar_index(:,:)
    logical,         allocatable :: detect_solvent_polar(:,:)
    integer                      :: natom, iatm
    integer                      :: atom_no
    character(6)                 :: residue_name, atomcls
    integer                      :: pre_resno
    integer                      :: n_partner
    type(s_partner), allocatable :: idx_to_atom(:)
    real(wp)                     :: mass


    n_solvent = size(option%solvent_list)
    allocate(detect_solvent(n_solvent), n_polar(n_solvent))

    ! count polar atoms in solvent molecules
    !
    n_polar(:) = 0
    detect_solvent(:) = .false.
    pre_resno = -99999

    natom = molecule%num_atoms
    do iatm = 1, natom
      do i_solvent = 1, n_solvent
        if (pre_resno /= molecule%residue_no(iatm) .and.  &
            n_polar(i_solvent) > 0) then

          detect_solvent(i_solvent) = .true.
        end if

        if (molecule%residue_name(iatm) == option%solvent_list(i_solvent) .and. &
            .not. detect_solvent(i_solvent)) then

          mass    = molecule%mass(iatm)
          atomcls = molecule%atom_cls_name(iatm)
          if (is_polar_atom(mass, atomcls)) then
            n_polar(i_solvent) = n_polar(i_solvent) + 1
          end if
        end if
      end do

      if (count(detect_solvent) == n_solvent)  &
        exit

      pre_resno = molecule%residue_no(iatm)
    end do

    ! setup
    !
    allocate(detect_solvent_polar(maxval(n_polar), n_solvent),  &
             solvent_polar_index(maxval(n_polar), n_solvent))

    pre_resno = -99999
    solvent_polar_index(:,:) = 0
    detect_solvent_polar(:,:) = .false.

    ! extract polar atoms only from the selected atom group
    !
    n_partner = 0
    natom = size(polar_atom)
    allocate(idx_to_atom(natom), partner_idx(natom))

    do iatm = 1, natom
      atom_no      = polar_atom(iatm)%atom_no
      atomcls      = molecule%atom_cls_name(atom_no)
      residue_name = molecule%residue_name(atom_no)
      mass         = molecule%mass(atom_no)

      if (.not. is_polar_atom(mass, atomcls))  &
        call error_msg('Setup_HB_partner_List> ERROR: non-polar atom is included')

      if (molecule%residue_no(atom_no) /= pre_resno) then
        i_polar = 0
      end if

      ! solvent atoms
      !
      do i_solvent = 1, n_solvent
        if (residue_name == option%solvent_list(i_solvent)) then
          i_polar = i_polar + 1

          if (.not. detect_solvent_polar(i_polar, i_solvent)) then
            n_partner = n_partner + 1
            solvent_polar_index(i_polar, i_solvent) = n_partner

            partner_idx(iatm) = n_partner
            detect_solvent_polar(i_polar, i_solvent) = .true.

            idx_to_atom(n_partner)%solvent = .true.
            idx_to_atom(n_partner)%atom_no = atom_no

          else
            partner_idx(iatm) = solvent_polar_index(i_polar, i_solvent)

          end if

          goto 10
        end if
      end do

      ! non-solvent atoms
      !
      n_partner = n_partner + 1
      partner_idx(iatm) = n_partner
      idx_to_atom(n_partner)%solvent = .false.
      idx_to_atom(n_partner)%atom_no = atom_no

 10   continue

      pre_resno = molecule%residue_no(atom_no)

    end do

    if (n_partner > 0) then
      allocate(partner_atom(n_partner))
      do iatm = 1, n_partner
        partner_atom(iatm)%solvent = idx_to_atom(iatm)%solvent
        partner_atom(iatm)%atom_no = idx_to_atom(iatm)%atom_no
      end do
    end if

    ! deallocate
    !
    deallocate(idx_to_atom)
    deallocate(detect_solvent, n_polar)
    deallocate(detect_solvent_polar, solvent_polar_index)

#ifdef DEBUG
    do iatm = 1, n_partner
      atom_no = partner_atom(iatm)%atom_no
      write(MsgOut,'(2I7,2(1x,A6))') iatm, atom_no,  &
                                     molecule%atom_name(atom_no),  &
                                     molecule%residue_name(atom_no)
    end do
#endif

    return

  end subroutine setup_hb_partner_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    examine_hbond
  !> @brief        examine the two polar atoms H-bond
  !! @authors      DM
  !! @param[in]    analysis_atom : atom you want to analyze
  !! @param[in]    target atom   : target atom for the H-bond analysis
  !! @param[in]    trajectory : trajectory information
  !! @param[in]    option     : option information
  !! @param[out]   hbond      : the two atom is H-bonded or not
  !! @param[inout] hb_list    : H-bond partner information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine examine_hbond(analysis_atom, target_atom, trajectory, option, &
                           hbond, hb_list)

    ! formal arguments
    type(s_polar_atom),         intent(in)    :: analysis_atom
    type(s_polar_atom),         intent(in)    :: target_atom
    type(s_trajectory), target, intent(in)    :: trajectory
    type(s_option),             intent(in)    :: option
    logical,                    intent(out)   :: hbond
    type(s_hb_info),            intent(inout) :: hb_list

    ! local variables
    real(wp), pointer :: coord(:,:), pbc_box(:,:)
    integer           :: ic
    integer           :: iatm, jatm, hatm
    integer           :: ih, ih_num, jh_num
    real(wp)          :: a_crd(3), t_crd(3), h_crd(3)
    real(wp)          :: diff(3), move(3)
    real(wp)          :: dist, dha_angle, hda_angle
    logical           :: pbc_correct


    pbc_correct = (option%boundary_type == BoundaryTypePBC)

    hbond   = .false.
    coord   => trajectory%coord
    pbc_box => trajectory%pbc_box

    ! atom number of the analysis and the target atoms
    iatm   = analysis_atom%atom_no
    jatm   =   target_atom%atom_no

    ! the numbers of hydrogen atoms bonded to the anlysis/target atoms
    ih_num = analysis_atom%num_hydrogen
    jh_num =   target_atom%num_hydrogen

    ! both polar atoms are H-bond acceptor
    if (ih_num == 0 .and. jh_num == 0)  &
      return

    ! exclude the case that iatm and jatm stand for the same atom
    if (iatm == jatm) &
      return

    ! move the target atom into the neighboring cell
    do ic = 1,3
      if (pbc_correct) then
        a_crd(ic) = coord(ic, iatm)  &
                    - floor(coord(ic, iatm) / pbc_box(ic, ic)) * pbc_box(ic, ic)
        t_crd(ic) = coord(ic, jatm)  &
                    - floor(coord(ic, jatm) / pbc_box(ic, ic)) * pbc_box(ic, ic)

        diff(ic) = a_crd(ic) - t_crd(ic)
        diff(ic) = diff(ic) - anint(diff(ic) / pbc_box(ic, ic)) * pbc_box(ic, ic)

      else
        a_crd(ic) = coord(ic, iatm)
        t_crd(ic) = coord(ic, jatm)

        diff(ic) = a_crd(ic) - t_crd(ic)
      end if

      if (abs(diff(ic)) > option%hb_distance)  &
        return
    end do

    ! judge by the H-bond distance
    dist = sqrt(dot_product(diff, diff))
    if (dist > option%hb_distance) &
      return


    ! compute the D-H .. A angle
    ! the case that iatm (analysis_atom) is the H-bond donor
    !
    if (pbc_correct) then
      do ic = 1, 3
        diff(ic) = a_crd(ic) - t_crd(ic)

        if (abs(diff(ic)) > (pbc_box(ic, ic) * half)) then
          move(ic) = sign(pbc_box(ic, ic), diff(ic))
          t_crd(ic) = t_crd(ic) + move(ic)
        end if

      end do
    end if

    do ih = 1, analysis_atom%num_hydrogen
      hatm = analysis_atom%hydrogen_atom(ih)

      if (pbc_correct) then
        do ic = 1, 3
          move(ic)  = a_crd(ic) - coord(ic, iatm)
          h_crd(ic) = coord(ic, hatm) + move(ic)
        end do
      else
        h_crd(:) = coord(:, hatm)
      end if

      dha_angle = compute_ang(a_crd, h_crd, t_crd)
      hda_angle = compute_ang(h_crd, a_crd, t_crd)

      if (dha_angle > option%dha_angle .and.  &
          hda_angle < option%hda_angle) then

        hb_list%partner   = jatm
        hb_list%hb_dist   = dist
        hb_list%dha_angle = dha_angle
        hb_list%hda_angle = hda_angle

        hbond = .true.
        exit
       end if
    end do

    ! the case that jatm (target_atom) is the H-bond donor
    if (pbc_correct) then
      do ic = 1,3
        diff(ic) = t_crd(ic) - a_crd(ic)

        if (abs(diff(ic)) > (pbc_box(ic, ic) * half)) then
          move(ic) = sign(pbc_box(ic, ic), diff(ic))
          a_crd(ic) = a_crd(ic) + move(ic)
        end if
      end do
    end if

    do ih = 1, target_atom%num_hydrogen
      hatm = target_atom%hydrogen_atom(ih)

      if (pbc_correct) then
        do ic = 1, 3
          move(ic)  = t_crd(ic) - coord(ic, jatm)
          h_crd(ic) = coord(ic, hatm) + move(ic)
        end do
      else
        h_crd(:) = coord(:, hatm)
      end if

      dha_angle = compute_ang(t_crd, h_crd, a_crd)
      hda_angle = compute_ang(h_crd, t_crd, a_crd)

      if (dha_angle > option%dha_angle .and.  &
          hda_angle < option%hda_angle) then

        if (hbond) then
          hbond = .false.
          write(MsgOut,'(i0," and ",i0," is in the geometry of D-H .. H-D.")') iatm, jatm
          exit
        end if

        hb_list%partner   = jatm
        hb_list%hb_dist   = dist
        hb_list%dha_angle = dha_angle
        hb_list%hda_angle = hda_angle

        hbond = .true.
        exit
      end if
    end do

    coord   => null()
    pbc_box => null()

    return

  end subroutine examine_hbond

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      is_polar_atom
  !> @brief        judge the atom is a polar atom or not
  !! @authors      DM
  !! @return       is_polar_atom
  !! @param[in]    atom_type
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
              
  function is_polar_atom(mass, atom_type)

    ! formal variables
    real(wp),     intent(in)    :: mass
    character(6), intent(inout) :: atom_type
    logical                     :: is_polar_atom


    if (atomic_number(mass, atom_type) == 7 .or.  &
        atomic_number(mass, atom_type) == 8) then

      is_polar_atom = .true.

    else
      is_polar_atom = .false.

    end if

    return

  end function is_polar_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    print_output_info
  !> @brief        print detailed output information
  !! @authors      TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8 

  subroutine print_output_info(output, option)

    ! formal arguments
    type(s_output),          intent(in) :: output
    type(s_option),          intent(in) :: option


    write(MsgOut,'(A)') ''
    write(MsgOut,'(A)') 'Analyze> Detailed information in the output files'
    write(MsgOut,'(A)') ''

    if (output%hb_listfile /= '') then
      write(MsgOut,'(A)') '  [hb_listfile] ' // trim(output%hb_listfile)
      write(MsgOut,'(A)') '    Column 1: snapshot index'
      write(MsgOut,'(A)') '    Column 2: "|"'
      write(MsgOut,'(A)') '    Column 3: atom name of the analysis atom'
      write(MsgOut,'(A)') '    Column 4: residue name of the analysis atom'
      write(MsgOut,'(A)') '    Column 5: residue number of the analysis atom'
      write(MsgOut,'(A)') '    Column 6: segment name of the analysis atom'
      write(MsgOut,'(A)') '    Column 7: ".."'
      write(MsgOut,'(A)') '    Column 8: atom name of the target atom'
      write(MsgOut,'(A)') '    Column 9: segment name of the target atom'
      write(MsgOut,'(A)') '    Column10: residue number of the target atom'
      write(MsgOut,'(A)') '    Column11: segment name of the target atom'
      write(MsgOut,'(A)') '    Column12: "|"'
      write(MsgOut,'(A)') '    Couumn13: H-bond distance (angstrom)'
      write(MsgOut,'(A)') '    Column14: H-bond H-D..A angle (degree)'
      write(MsgOut,'(A)') '    Column15: H-bond D-H..A angle (degree)'
      write(MsgOut,'(A)') ''
    end if

    select case (option%output_type)

    case (HBOutputModeCountSnap)

      write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
      write(MsgOut,'(A)') '    output type is "count_snap"'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '    Column 1: snapshot index'
      write(MsgOut,'(A)') '    Column 2: the number of H-bonds formed in the snapshot'
      write(MsgOut,'(A)') ''

    case (HBOutputModeCountAtom)
      write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
      write(MsgOut,'(A)') '    output type is "count_atom"'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '    Column 1: count of H-bonds formed between the two atoms in the input trajectory'
      write(MsgOut,'(A)') '    Column 2: "|"'
      write(MsgOut,'(A)') '    Column 3: atom name of the analysis atom'
      write(MsgOut,'(A)') '    Column 4: residue name of the analysis atom'
      write(MsgOut,'(A)') '    Column 5: residue number of the analysis atom'
      write(MsgOut,'(A)') '    Column 6: segment name of the analysis atom'
      write(MsgOut,'(A)') '    Column 7: ".."'
      write(MsgOut,'(A)') '    Column 8: atom name of the target atom'
      write(MsgOut,'(A)') '    Column 9: segment name of the target atom'
      write(MsgOut,'(A)') '    Column10: residue number of the target atom (not output for "solvent_list")'
      write(MsgOut,'(A)') '    Column11: segment name of the target atom (not output for "solvent_list")'
      write(MsgOut,'(A)') ''

    case (HBOutputModeLifetime)
      write(MsgOut,'(A)') '  [outfile] ' // trim(output%outfile)
      write(MsgOut,'(A)') '    output type is "lifetime"'
      write(MsgOut,'(A)') ''
      write(MsgOut,'(A)') '    Column 1: "snapshot"'
      write(MsgOut,'(A)') '    Column 2: the last snapshot number of the H-bond'
      write(MsgOut,'(A)') '    Column 3: ":"'
      write(MsgOut,'(A)') '    Column 4: lifetime of the continuously formed H-bond '
      write(MsgOut,'(A)') '              (the count of snapshots in which the H-bond was retained)'
      write(MsgOut,'(A)') '    Column 5: "|"'
      write(MsgOut,'(A)') '    Column 6: atom name of the analysis atom'
      write(MsgOut,'(A)') '    Column 7: residue name of the analysis atom'
      write(MsgOut,'(A)') '    Column 8: residue number of the analysis atom'
      write(MsgOut,'(A)') '    Column 9: segment name of the analysis atom'
      write(MsgOut,'(A)') '    Column10: ".."'
      write(MsgOut,'(A)') '    Column11: atom name of the target atom'
      write(MsgOut,'(A)') '    Column12: segment name of the target atom'
      write(MsgOut,'(A)') '    Column13: residue number of the target atom (not output for "solvent list")'
      write(MsgOut,'(A)') '    Column14: segment name of the target atom (not output for "solvent list")'
      write(MsgOut,'(A)') ''

    end select

    return

  end subroutine print_output_info

end module hb_analyze_mod
