!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   hbond_analyze_mod
!> @brief   run hbond analysis
!! @authors Daisuke Matsuoka (DM), Isseki Yu (IY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module hbond_analyze_mod

  use hbond_option_str_mod
  use sa_domain_mod
  use sa_boundary_mod
  use sa_tool_mod
  use sa_domain_str_mod
  use sa_boundary_str_mod
  use sa_ensemble_str_mod
  use sa_option_str_mod
  use trajectory_str_mod
  use output_str_mod
  use fileio_trj_mod
  use fileio_mod
  use molecules_str_mod
  use measure_mod
  use select_atoms_str_mod
  use atom_libs_mod
  use messages_mod
  use constants_mod
  use mpi_parallel_mod
  use string_mod
#ifdef DEBUG
  use molecules_mod
  use fileio_pdb_mod
#endif
#ifdef HAVE_MPI_GENESIS
  use mpi
#endif
#ifdef OMP
  use omp_lib
#endif

  implicit none
#ifdef HAVE_MPI_GENESIS
#ifdef MSMPI
!GCC$ ATTRIBUTES DLLIMPORT :: MPI_BOTTOM, MPI_IN_PLACE
#endif
#endif
  private

  ! structures
  type, private :: s_polar_info
    integer :: num_hydrogen
    integer :: hydrogen_atom(4)    ! e.g. NH4+ has four 'H-bond' donor.
  end type s_polar_info

  type, private :: s_pl_list
    integer,            allocatable :: idx(:)
    type(s_polar_info), allocatable :: info(:)
  end type s_pl_list

  type, private :: s_partner
    logical :: solvent
    integer :: atom_no
  end type s_partner

  type, private :: s_hb_info
    real(sp) :: hb_dist
    real(sp) :: dha_angle
    real(sp) :: hda_angle
  end type s_hb_info

  ! subroutines
  public  :: run_hbond
  private :: get_polar_atom
  private :: setup_hb_partner_list
  private :: examine_hbond
  private :: is_polar_atom
  private :: pio_check_ranked_file
  private :: pio_get_ranked_filename
  private :: check_options

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    analyze
  !> @brief        analyze trajectories to caluclate hbond distribution
  !! @authors      DM, IY
  !! @param[inout] molecule     : molecule information
  !! @param[in]    trj_list     : trajectory file list information
  !! @param[in]    output       : output information
  !! @param[inout] option       : option information
  !! @param[in]    hbond_option : hbond_option information
  !! @param[inout] trajectory   : trajectory information
  !! @param[inout] boundary     : boundary information
  !! @param[inout] domain       : domain information
  !! @param[in]    ensemble     : ensemble information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine run_hbond(molecule, trj_list, output, option,  &
                       hbond_option, trajectory, boundary, domain, ensemble)

    ! formal arguments
    type(s_molecule), target, intent(inout) :: molecule
    type(s_trj_list),         intent(in)    :: trj_list
    type(s_output),           intent(in)    :: output
    type(s_option),   target, intent(inout) :: option
    type(s_hbond_option),     intent(in)    :: hbond_option
    type(s_trajectory),       intent(inout) :: trajectory
    type(s_boundary),         intent(inout) :: boundary
    type(s_domain),           intent(inout) :: domain
    type(s_ensemble),         intent(in)    :: ensemble

    ! local variables
    type(s_trj_file)             :: trj_in
    real(wp)                     :: cofm(3), move(3)
    real(wp)                     :: box_size(3)
    integer                      :: ifile, istep, totstep, nfile, i, j
    integer                      :: out_unit, pio_unit
    integer                      :: alloc_stat, dealloc_stat
    integer                      :: my_id
    integer                      :: igroup, ngroup
    integer                      :: gid_analysis, gid_target
    integer                      :: natom_selec
    integer                      :: iatom, jatom, natom
    integer                      :: inside, icount
    integer                      :: alloc_step
    integer                      :: icell
    integer                      :: iatom_analysis, iatom_target
    integer                      :: solvent_idx
    integer                      :: hb_total_dom, hb_total
    integer                      :: a_atm, t_atm
    character(len=MaxFilename)   :: pio_filename

#ifdef DEBUG
    type(s_pdb)                  :: check_trj
#endif

    type(s_pl_list), allocatable, target :: polar_list(:)
    type(s_partner), allocatable :: partner_atom(:)
    type(s_parray),  allocatable :: selec_atom(:)
    type(s_hb_info), allocatable :: hb_info(:)
    logical,         allocatable :: hbond(:)
    integer,         allocatable :: partner_idx(:)
    integer,         allocatable :: offset(:), natom_group(:)
    integer,         allocatable :: analysis_in_domain(:), target_in_bound(:)
    integer,         allocatable :: hb_count(:,:)

    real(wp),            pointer :: atom_coord(:,:)
    integer,             pointer :: numr(:)
    character(len=4),    pointer :: nama(:)
    character(len=6),    pointer :: namr(:)
    character(len=4),    pointer :: seg(:)


    ! formats
100 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',1x,a4,1x,a6,1x,i7,1x,a4' | ',F6.3,2F9.3)
102 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',1x,a4,1x,a6)
104 format(i10,' | ',a4,1x,a6,1x,i7,1x,a4,' .. ',1x,a4,1x,a6,1x,i7,1x,a4)
106 format('snapshot',i10,' : ',i10,' | ',i7,1x,a4,1x,a6,1x,i7,1x,a4,' .. ',i7,1x,a4,1x,a6,1x,i7,1x,a4)

    ! check options
    !
    totstep      = 0
    call check_options(molecule, boundary, domain, trj_list, option, &
                       hbond_option, totstep)

    nama       => molecule%atom_name
    namr       => molecule%residue_name
    numr       => molecule%residue_no
    seg        => molecule%segment_name
    atom_coord => molecule%atom_coord

    natom        = molecule%num_atoms
    nfile        = size(trj_list%md_steps)
    ngroup       = size(option%analysis_atoms)
    gid_analysis = hbond_option%analysis_atom
    gid_target   = hbond_option%target_atom

    ! allocation
    !
    allocate(polar_list(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc

    allocate(selec_atom(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc

    allocate(natom_group(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc

    allocate(offset(ngroup), stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc


    ! list polar atoms in groups
    !
    offset(:)   = 0
    natom_selec = 0
    do igroup = 1, ngroup
      call get_polar_atom(molecule, option%analysis_atoms(igroup), &
                          polar_list(igroup))

      if (igroup == gid_target) then
        call setup_hb_partner_list(molecule, hbond_option, polar_list(igroup),&
                                   partner_idx, partner_atom)
      end if

      selec_atom(igroup)%idx => polar_list(igroup)%idx
      natom_group(igroup)    = size(polar_list(igroup)%idx)
      offset(igroup)         = natom_selec 
      natom_selec            = natom_selec + natom_group(igroup) 
    end do

    domain%num_group      = ngroup
    domain%num_atom_selec = natom_selec

    ! allocation and initialization
    !
    allocate(hb_count(size(partner_atom), natom_group(gid_analysis)), &
             stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc


    do i = 1, natom_group(gid_analysis)
      do j = 1, size(partner_atom)
        hb_count(j,i) = 0
      end do
    end do

    ! open ouput file
    !
    if (output%hb_listfile .ne. '') then

      if (pio_check_ranked_file(output%hb_listfile)) then

        pio_filename = pio_get_ranked_filename(output%hb_listfile)
        call open_file(pio_unit, pio_filename, IOFileOutputNew)

        write(pio_unit,'(A)') '# This file was created with hbond_analysis'
        write(pio_unit,'(A)') '# '
        write(pio_unit,'(A)') '# Column 1: snapshot index'
        write(pio_unit,'(A)') '# Column 2: "|"'
        write(pio_unit,'(A)') '# Column 3: atom name of the analysis atom'
        write(pio_unit,'(A)') '# Column 4: residue name of the analysis atom'
        write(pio_unit,'(A)') '# Column 5: residue number of the analysis atom'
        write(pio_unit,'(A)') '# Column 6: segment name of the analysis atom'
        write(pio_unit,'(A)') '# Column 7: ".."'
        write(pio_unit,'(A)') '# Column 8: atom name of the target atom'
        write(pio_unit,'(A)') '# Column 9: segment name of the target atom'
        write(pio_unit,'(A)') '# Column10: residue number of the target atom'
        write(pio_unit,'(A)') '# Column11: segment name of the target atom'
        write(pio_unit,'(A)') '# Column12: "|"'
        write(pio_unit,'(A)') '# Couumn13: H-bond distance (angstrom)'
        write(pio_unit,'(A)') '# Column14: H-bond H-D..A angle (degree)'
        write(pio_unit,'(A)') '# Column15: H-bond D-H..A angle (degree)'
        write(pio_unit,'(A)') ''

      else
        call error_msg('pio_output is necessary if you want to output H-Bond list per snapshot.')
      end if

    end if

    if (output%txtfile .ne. '') then
      if (main_rank) then
        call open_file(out_unit, output%txtfile, IOFileOutputNew)

        select case (hbond_option%output_type)

        case (HBOutputModeCountSnap)
          write(out_unit,'(A)') '# This file was created with hbond_analysis'
          write(out_unit,'(A)') '# output type is "count_snap"'
          write(out_unit,'(A)') '# '
          write(out_unit,'(A)') '# Column 1: snapshot index'
          write(out_unit,'(A)') '# Column 2: the number of H-bonds formed in the snapshot'
          write(out_unit,'(A)') ''

        case (HBOutputModeCountAtom)
          write(out_unit,'(A)') '# This file was created with hbond_analysis'
          write(out_unit,'(A)') '# output type is "count_atom"'
          write(out_unit,'(A)') '#'
          write(out_unit,'(A)') '# Column 1: count of H-bonds formed between the two atoms in the trajectory'
          write(out_unit,'(A)') '# Column 2: "|"'
          write(out_unit,'(A)') '# Column 3: atom name of the analysis atom'
          write(out_unit,'(A)') '# Column 4: residue name of the analysis atom'
          write(out_unit,'(A)') '# Column 5: residue number of the analysis atom'
          write(out_unit,'(A)') '# Column 6: segment name of the analysis atom'
          write(out_unit,'(A)') '# Column 7: ".."'
          write(out_unit,'(A)') '# Column 8: atom name of the target atom'
          write(out_unit,'(A)') '# Column 9: segment name of the target atom'
          write(out_unit,'(A)') '# Column10: residue number of the target atom (not output for "solvent_list")'
          write(out_unit,'(A)') '# Column11: segment name of the target atom (not output for "solvent_list")'
          write(out_unit,'(A)') ''

        end select

      end if
    end if

    ! analysis loop
    !
    do ifile = 1, nfile

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

          totstep = totstep + 1
          if (main_rank) then
            write(MsgOut,*) '      number of structures = ', totstep
          end if

          do iatom = 1, molecule%num_atoms
            molecule%atom_coord(1:3,iatom) = trajectory%coord(1:3, iatom)
          end do

          box_size(1) = trajectory%pbc_box(1,1)
          box_size(2) = trajectory%pbc_box(2,2)
          box_size(3) = trajectory%pbc_box(3,3)

          ! centering
          !
          if (hbond_option%recenter /= 0) then
            call get_cofm(molecule, selec_atom, hbond_option%recenter, cofm)
            move(1:3) = -cofm(1:3)
            call shift_molecule(move, molecule)
          end if

          ! wrap molecule
          !
          call wrap_molecule(option, trajectory, ensemble, &
                             boundary, molecule)

          ! setup domains
          !
          call setup_atom(option, ensemble, &
                          trajectory, selec_atom, molecule, boundary, domain)

#ifdef DEBUG
          if (main_rank) then
            call export_molecules(molecule, pdb = check_trj)
            call output_pdb('hbond_check_trj.pdb', check_trj)
          end if
#endif

          ! determinie analysis atoms in domain
          !
          do alloc_step = 1, 2
            icount = 0

            do iatom = 1, natom_group(gid_analysis)
              icell = domain%selec_atom2cell(iatom + offset(gid_analysis))

              if (icell > 0) then
                inside = domain%cell_g2l(icell)

                if (inside > 0) then
                  icount = icount + 1

                  if (alloc_step == 2) then
                    analysis_in_domain(icount) = iatom
                  end if
                end if
              end if
            end do

            if (alloc_step == 1) then
              if (allocated(analysis_in_domain)) then
                deallocate(analysis_in_domain, stat = dealloc_stat)
                if (dealloc_stat /= 0) &
                  call error_msg_dealloc

              end if

              allocate(analysis_in_domain(icount), stat = alloc_stat)
              if (alloc_stat /= 0) &
                call error_msg_alloc

              analysis_in_domain(:) = 0
            end if
          end do

#ifdef DEBUG
          write(MsgOut,*) 'analysis atoms', &
               my_world_rank, size(analysis_in_domain), &
               size(selec_atom(gid_analysis)%idx)
#endif

          ! determinie target atoms in domain + boundary
          !
          do alloc_step = 1, 2
            icount = 0

            do iatom = 1, natom_group(gid_target)
              icell = domain%selec_atom2cell(iatom + offset(gid_target))

              if (icell > 0) then
                inside = domain%cell_g2l(icell) + domain%cell_g2b(icell)

                if (inside > 0) then
                  icount = icount + 1

                  if (alloc_step == 2) then
                    target_in_bound(icount) = iatom
                  end if
                end if
              end if
            end do

            if (alloc_step == 1) then
              if (allocated(target_in_bound)) then
                deallocate(target_in_bound, &
                           hb_info, &
                           hbond,   &
                           stat = dealloc_stat)
                if (dealloc_stat /= 0) &
                  call error_msg_dealloc

              end if

              allocate(target_in_bound(icount), &
                       hb_info(icount), &
                       hbond(icount), stat = alloc_stat)
              if (alloc_stat /= 0) &
                call error_msg_alloc

              target_in_bound(:) = 0
            end if
          end do

#ifdef DEBUG
          write(MsgOut,*) 'target atoms', &
               my_world_rank, size(target_in_bound),  &
               size(selec_atom(gid_target)%idx)
#endif

          hb_total_dom = 0

          do iatom = 1, size(analysis_in_domain)

            iatom_analysis = analysis_in_domain(iatom)

            !$omp parallel default(shared)              &
            !$omp private (my_id, jatom, solvent_idx)   &
            !$omp private (t_atm, iatom_target)         &
            !$omp reduction(+ : hb_total_dom)
#ifdef OMP
            my_id = omp_get_thread_num()
#else
            my_id = 0
#endif

            do jatom = my_id+1, size(target_in_bound), nthread

              iatom_target   = target_in_bound(jatom)

              call examine_hbond(polar_list(gid_analysis)%idx(iatom_analysis), &
                                 polar_list(gid_target)%idx(iatom_target),     &
                                 polar_list(gid_analysis)%info(iatom_analysis),&
                                 polar_list(gid_target)%info(iatom_target),    &
                                 atom_coord, box_size, hbond_option,           &
                                 hbond(jatom), hb_info(jatom))

              select case (hbond_option%output_type)

                case (HBOutputModeCountSnap)

                  if (hbond(jatom)) then
                    hb_total_dom = hb_total_dom + 1
                  end if

                case (HBOutputModeCountAtom)

                  if (hbond(jatom)) then
                    solvent_idx = partner_idx(iatom_target)
                    !$omp critical
                    hb_count(solvent_idx, iatom_analysis)   &
                            = hb_count(solvent_idx, iatom_analysis) + 1
                    !$omp end critical
                  end if

              end select

            end do
            !$omp end parallel

            do jatom = 1, size(target_in_bound)

              iatom_target = target_in_bound(jatom)

              if (pio_check_ranked_file(output%hb_listfile)) then
                if (hbond(jatom)) then
                  a_atm = polar_list(gid_analysis)%idx(iatom_analysis)
                  t_atm = polar_list(gid_target)%idx(iatom_target)

                  write(pio_unit, 100) &
                       totstep, &
                       nama(a_atm), namr(a_atm), numr(a_atm), seg(a_atm), &
                       nama(t_atm), namr(t_atm), numr(t_atm), seg(t_atm), &
                       hb_info(jatom)%hb_dist,   &
                       hb_info(jatom)%dha_angle, &
                       hb_info(jatom)%hda_angle
                end if
              end if

            end do

          end do

          ! output results
          if ((hbond_option%output_type == HBOutputModeCountSnap) .and. &
              (output%txtfile .ne. '')) then

#ifdef HAVE_MPI_GENESIS
            call mpi_reduce(hb_total_dom, hb_total, 1, mpi_integer, &
                            mpi_sum, 0, mpi_comm_country, ierror)
#else
            hb_total = hb_total_dom
#endif
            if (main_rank) then
              write(out_unit, '(I10,1x,I10)') totstep, hb_total
            end if

          end if

        end if
      end do

      ! close trajectory file
      !
      call close_trj(trj_in)

    end do

    ! output results
    if ((hbond_option%output_type == HBOutputModeCountAtom) .and.  &
        (output%txtfile .ne. '')) then

#ifdef DEBUG
      write(6,*) 'rank = ', my_country_rank, sizeof(hb_count)/4, &
                 size(partner_atom), natom_group(gid_analysis)
#endif

#ifdef HAVE_MPI_GENESIS
      if (main_rank) then
        do i = 1, natom_group(gid_analysis)
          call mpi_reduce(mpi_in_place, hb_count(:,i), size(hb_count(:,i)), &
                          mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        end do
      else
        do i = 1, natom_group(gid_analysis)
          call mpi_reduce(hb_count(:,i), 0, size(hb_count(:,i)), &
                          mpi_integer, mpi_sum, 0, mpi_comm_country, ierror)
        end do
      end if

#else

#ifdef DEBUG
        do i = 1, natom_group(gid_analysis)
          write(0,*) (hb_count(j, i) , j = 1, size(partner_atom))
        end do
#endif

#endif

      if (main_rank) then
        do iatom = 1, natom_group(gid_analysis)
          do jatom = 1, size(partner_atom)
            if (hb_count(jatom, iatom) > 0) then
              a_atm = polar_list(gid_analysis)%idx(iatom)
              t_atm = partner_atom(jatom)%atom_no

              if (partner_atom(jatom)%solvent) then
                write(out_unit, 102) &
                     hb_count(jatom, iatom), &
                     nama(a_atm), namr(a_atm), numr(a_atm), seg(a_atm), &
                     nama(t_atm), namr(t_atm)
              else
                write(out_unit, 104) &
                     hb_count(jatom, iatom), &
                     nama(a_atm), namr(a_atm), numr(a_atm), seg(a_atm), &
                     nama(t_atm), namr(t_atm), numr(t_atm), seg(t_atm)
              end if
            end if
          end do
        end do
      end if

    end if

    if (main_rank) then
      write(MsgOut,*) ' '
    end if

    ! close output file
    !
    if (output%hb_listfile .ne. '')  &
      call close_file(pio_unit)

    if (main_rank) then
      call close_file(out_unit)
    end if

    ! deallocation
    !
    deallocate(selec_atom,  stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(polar_list,  stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(natom_group, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(offset,      stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(partner_idx, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(partner_atom, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(analysis_in_domain, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(target_in_bound, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(hb_count, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(hb_info, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    deallocate(hbond, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    nullify(atom_coord, nama, namr, numr, seg)

    return

  end subroutine run_hbond

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
              
  subroutine get_polar_atom(molecule, select_atom, polar_list)

    ! formal variables
    type(s_molecule),        intent(in)    :: molecule
    type(s_selatoms),        intent(in)    :: select_atom
    type(s_pl_list),         intent(out)   :: polar_list

    ! local variables
    integer                  :: alloc_stat, dealloc_stat
    integer                  :: npolar, ipolar
    integer                  :: iatm, natom, atom_no, ibond
    integer                  :: bond_atom1, bond_atom2
    character(len=6)         :: atomcls1, atomcls2 
    character(len=3)         :: atom_element

    integer, allocatable     :: atom_to_idx(:)


    ! extract polar atoms only from the selected atom group
    !
    npolar = 0
    natom  = size(select_atom%idx)

    do iatm = 1, natom
      atom_no      = select_atom%idx(iatm)
      atom_element = molecule%atom_cls_name(atom_no)(1:3)

      if (is_polar_atom(atom_element)) then
        npolar = npolar + 1
      end if

    end do

    if (allocated(polar_list%idx)) then
      deallocate(polar_list%idx, stat = dealloc_stat)
      if (dealloc_stat /= 0)  call error_msg_dealloc
    end if

    if (allocated(polar_list%info)) then
      deallocate(polar_list%info, stat = dealloc_stat)
      if (dealloc_stat /= 0)  call error_msg_dealloc
    end if

    allocate(polar_list%idx(npolar),  &
             polar_list%info(npolar), &
             stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc

    polar_list%idx(:)               = 0
    polar_list%info(:)%num_hydrogen = 0

    allocate(atom_to_idx(molecule%num_atoms), stat = alloc_stat)
    if (alloc_stat /= 0)  call error_msg_alloc
    atom_to_idx(:) = -99999

    npolar = 0
    do iatm = 1, natom
      atom_no      = select_atom%idx(iatm)
      atom_element = molecule%atom_cls_name(atom_no)(1:3)

      if (is_polar_atom(atom_element)) then

        npolar = npolar + 1
        polar_list%idx(npolar) = atom_no
        atom_to_idx(atom_no) = npolar

      end if
    end do

    ! get hydrogen atoms bonded to the polar atoms
    !
    do ibond = 1, molecule%num_bonds
      bond_atom1 = molecule%bond_list(1, ibond)
        atomcls1 = molecule%atom_cls_name(bond_atom1)
      
      bond_atom2 = molecule%bond_list(2, ibond)
        atomcls2 = molecule%atom_cls_name(bond_atom2)

      if (atomic_number_by_name(atomcls1(1:3)) /= 1 .and. &
          atomic_number_by_name(atomcls2(1:3)) /= 1)      &
        cycle

      if (atom_to_idx(bond_atom1) > 0) then
        ipolar = atom_to_idx(bond_atom1)
        polar_list%info(ipolar)%num_hydrogen = &
             polar_list%info(ipolar)%num_hydrogen + 1

        iatm = polar_list%info(ipolar)%num_hydrogen
        polar_list%info(ipolar)%hydrogen_atom(iatm) = bond_atom2

      else if (atom_to_idx(bond_atom2) > 0) then
        ipolar = atom_to_idx(bond_atom2)
        polar_list%info(ipolar)%num_hydrogen = &
             polar_list%info(ipolar)%num_hydrogen + 1

        iatm = polar_list%info(ipolar)%num_hydrogen
        polar_list%info(ipolar)%hydrogen_atom(iatm) = bond_atom1

      end if

    end do

    deallocate(atom_to_idx, stat = dealloc_stat)
    if (dealloc_stat /= 0)  call error_msg_dealloc

    return

  end subroutine get_polar_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_hb_partner_list
  !> @brief        make H-bond partner list
  !! @authors      DM
  !! @param[in]    molecule     : molecule information
  !! @param[in]    hbond_option : option information
  !! @param[in]    polar_list   : H-bondable atom information
  !! @param[out]   parnter_idx  : polar atom index to partner index
  !! @param[out]   partner_atom : H-bond partner information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_hb_partner_list(molecule, hbond_option, polar_list, &
                                   partner_idx, partner_atom)

    ! formal arguments
    type(s_molecule),             intent(in)    :: molecule
    type(s_hbond_option),         intent(in)    :: hbond_option
    type(s_pl_list),              intent(in)    :: polar_list
    integer,         allocatable, intent(out)   :: partner_idx(:)
    type(s_partner), allocatable, intent(out)   :: partner_atom(:)

    ! local variable
    integer                       :: n_solvent, i_solvent
    integer                       :: i_polar
    integer                       :: natom, iatm
    integer                       :: atom_no
    integer                       :: pre_resno
    integer                       :: n_partner
    character(len=6)              :: residue_name
    character(len=3)              :: atom_element

    integer,         allocatable  :: npolar(:)
    logical,         allocatable  :: detect_solvent(:)
    integer,         allocatable  :: solvent_polar_index(:,:)
    logical,         allocatable  :: detect_solvent_polar(:,:)
    type(s_partner), allocatable  :: idx_to_atom(:)


    n_solvent = size(hbond_option%solvent_list)
    allocate(detect_solvent(n_solvent), npolar(n_solvent))

    ! count polar atoms in solvent molecules
    !
    npolar(:) = 0
    detect_solvent(:) = .false.
    pre_resno = -9999999

    natom = molecule%num_atoms
    do iatm = 1, natom
      do i_solvent = 1, n_solvent
        if (pre_resno /= molecule%residue_no(iatm) .and.  &
            npolar(i_solvent) > 0) then

          detect_solvent(i_solvent) = .true.
        end if

        if (molecule%residue_name(iatm) == hbond_option%solvent_list(i_solvent) .and. &
            .not. detect_solvent(i_solvent)) then

          if (is_polar_atom(molecule%atom_cls_name(iatm))) then
            npolar(i_solvent) = npolar(i_solvent) + 1
          end if
        end if
      end do

      if (count(detect_solvent) == n_solvent)  &
        exit

      pre_resno = molecule%residue_no(iatm)
    end do

    ! setup
    !
    allocate(detect_solvent_polar(n_solvent, maxval(npolar)),  &
             solvent_polar_index (n_solvent, maxval(npolar)))

    pre_resno = -99999
    solvent_polar_index(:,:) = 0
    detect_solvent_polar(:,:) = .false.

    ! extract polar atoms only from the selected atom group
    !
    n_partner = 0
    natom = size(polar_list%idx)
    allocate(idx_to_atom(natom), partner_idx(natom))

    do iatm = 1, natom
      atom_no      = polar_list%idx(iatm)
      atom_element = molecule%atom_cls_name(atom_no)(1:3)
      residue_name = molecule%residue_name(atom_no)

      if (.not. is_polar_atom(atom_element))  &
        call error_msg( &
        'Setup_HB_partner_List> ERROR: non-polar atom is included')

      if (molecule%residue_no(atom_no) /= pre_resno) then
        i_polar = 0
      end if

      ! solvent atoms
      !
      do i_solvent = 1, n_solvent
        if (residue_name == hbond_option%solvent_list(i_solvent)) then
          i_polar = i_polar + 1

          if (.not. detect_solvent_polar(i_solvent, i_polar)) then
            n_partner = n_partner + 1
            solvent_polar_index(i_solvent, i_polar) = n_partner

            partner_idx(iatm) = n_partner
            detect_solvent_polar(i_solvent, i_polar) = .true.

            idx_to_atom(n_partner)%solvent = .true.
            idx_to_atom(n_partner)%atom_no = atom_no

          else
            partner_idx(iatm) = solvent_polar_index(i_solvent, i_polar)

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

10    continue

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
    deallocate(detect_solvent, npolar)
    deallocate(detect_solvent_polar, solvent_polar_index)

    return

  end subroutine setup_hb_partner_list

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    examine_hbond
  !> @brief        examine the two polar atoms H-bond
  !! @authors      DM
  !! @param[in]    analysis_atom : atom you want to analyze
  !! @param[in]    target atom   : target atom for the H-bond analysis
  !! @param[in]    analysis_info : hydrogen atom information of analysis atoms
  !! @param[in]    target_info   : hydrogen atom information of target atoms
  !! @param[in]    atom_coord    : atomic coordinate information
  !! @param[in]    box_size      : box dimension information
  !! @param[in]    hbond_option  : H-bond option information
  !! @param[out]   hbond         : the two atom is H-bonded or not
  !! @param[inout] hb_list       : H-bond partner information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine examine_hbond(analysis_atom, target_atom,  &
                           analysis_info, target_info,  &
                           atom_coord, box_size, hbond_option, hbond, hb_list)

    ! formal arguments
    integer,              intent(in)    :: analysis_atom
    integer,              intent(in)    :: target_atom
    type(s_polar_info),   intent(in)    :: analysis_info
    type(s_polar_info),   intent(in)    :: target_info
    real(wp),             intent(in)    :: atom_coord(:,:)
    real(wp),             intent(in)    :: box_size(3)
    type(s_hbond_option), intent(in)    :: hbond_option
    logical,              intent(out)   :: hbond
    type(s_hb_info),      intent(inout) :: hb_list

    ! local variables
    real(wp)              :: a_crd(3), t_crd(3), h_crd(3)
    real(wp)              :: diff(3), move(3)
    real(wp)              :: dist, dha_angle, hda_angle
    integer               :: hatm
    integer               :: ih, ih_num, jh_num
    integer               :: icomp


    hbond = .false.

    ! the numbers of hydrogen atoms bonded to the anlysis/target atoms
    ih_num = analysis_info%num_hydrogen
    jh_num =   target_info%num_hydrogen

    ! both polar atoms are H-bond acceptor
    if (ih_num == 0 .and. jh_num == 0)  &
      return

    ! exclude the case that iatm and jatm stand for the same atom
    if (analysis_atom == target_atom) &
      return

    ! coordinates of the analysis and the target atoms
    a_crd(:) = atom_coord(:, analysis_atom)
    t_crd(:) = atom_coord(:,   target_atom)

    ! move the target atom into the neighboring cell
    move(:) = 0.0_wp
    diff(:) = a_crd(:) - t_crd(:)

    do icomp = 1, 3
      if (abs(diff(icomp)) > (box_size(icomp) * half)) then
        move(:) = sign(box_size(icomp), diff(icomp))
        t_crd(icomp) = t_crd(icomp) + move(icomp)
      end if
    end do

    ! judge by the H-bond distance
    diff(:) = a_crd(:) - t_crd(:)

    if (abs(diff(1)) > hbond_option%hb_distance  .or.  &
        abs(diff(2)) > hbond_option%hb_distance  .or.  &
        abs(diff(3)) > hbond_option%hb_distance)       &
      return
      
    dist = sqrt(dot_product(diff, diff))
    if (dist > hbond_option%hb_distance) &
      return

    ! compute the D-H .. A angle
    ! the case that iatm (analysis_atom) is the H-bond donor
    !
    do ih = 1, analysis_info%num_hydrogen
      hatm = analysis_info%hydrogen_atom(ih)
      h_crd(:) = atom_coord(:,hatm)

      dha_angle = compute_ang(a_crd, h_crd, t_crd)
      hda_angle = compute_ang(h_crd, a_crd, t_crd)

      if (dha_angle > hbond_option%dha_angle .and.  &
          hda_angle < hbond_option%hda_angle) then

        hb_list%hb_dist   = real(dist, sp)
        hb_list%dha_angle = real(dha_angle, sp)
        hb_list%hda_angle = real(hda_angle, sp)

        hbond = .true.
        exit
      end if
    end do

    ! the case that jatm (target_atom) is the H-bond donor
    do ih = 1, target_info%num_hydrogen
      hatm = target_info%hydrogen_atom(ih)
      h_crd(:) = atom_coord(:,hatm) + move(:)

      dha_angle = compute_ang(t_crd, h_crd, a_crd)
      hda_angle = compute_ang(h_crd, t_crd, a_crd)

      if (dha_angle > hbond_option%dha_angle .and.  &
          hda_angle < hbond_option%hda_angle) then

        if (hbond) then
          hbond = .false.
          write(MsgOut,'(i0," and ",i0," is in the geometry of D-H .. H-D.")') &
            analysis_atom, target_atom
          return
        end if

        hb_list%hb_dist   = real(dist, sp)
        hb_list%dha_angle = real(dha_angle, sp)
        hb_list%hda_angle = real(hda_angle, sp)

        hbond = .true.
        exit
      end if
    end do

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
              
  function is_polar_atom(atom_type)

    ! function
    logical                  :: is_polar_atom

    ! formal variables
    character(*),            intent(in) :: atom_type


    if (atomic_number_by_name(atom_type) == 7 .or.  &
        atomic_number_by_name(atom_type) == 8) then

      is_polar_atom = .true.

    else
      is_polar_atom = .false.

    end if

    return

  end function is_polar_atom

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_check_ranked_file
  !> @brief        check ranked file name or not
  !! @authors      NT
  !! @param[in]    filename : file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_check_ranked_file(filename)

    ! function
    logical                  :: pio_check_ranked_file

    ! formal arguments
    character(100),          intent(in) :: filename

    ! local variables
    integer                  :: ci1, ci2


    ! check filename
    !
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    pio_check_ranked_file = (ci1 /= 0 .and. ci2 /= 0 .and. ci1 < ci2)

    return

  end function pio_check_ranked_file

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Function      pio_get_ranked_filename
  !> @brief        get ranked file name
  !! @authors      NT
  !! @param[in]    filename : file name
  !! @param[out]   nplace   : number of places
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function pio_get_ranked_filename(filename, nplace)

    ! function
    character(100)           :: pio_get_ranked_filename

    ! formal arguments
    character(*),            intent(in) :: filename
    integer,       optional, intent(in) :: nplace

    ! local variables
    integer                  :: ci1, ci2, ip
    character(50)            :: fmt_str

    ! constants
    integer,       parameter :: Places (7) = &
         (/10, 100, 1000, 10000, 100000, 1000000, 10000000/)


    if (present(nplace)) then
      ip = nplace
    else
      do ip = 1, size(Places)
        if (nproc_world < Places(ip)) &
          exit
      end do
    end if

    ! check filename
    !
    pio_get_ranked_filename = filename
    ci1 = scan(filename, '(', .true.)
    ci2 = scan(filename, ')', .true.)

    if (ci1 == 0 .or. ci2 ==0 .or. ci1 > ci2) then
      call error_msg( &
      'Pio_Get_Ranked_Filename> Filename is not correctly ranked in [OUTPUT]')
    end if

    write(fmt_str,'(A,I0,A,I0,A)') '(A,I',ip,'.',ip,',A)'

    write(pio_get_ranked_filename,fmt=fmt_str)   &
         trim(filename(1:ci1-1)), my_world_rank, &
         trim(filename(ci2+1:))

    return

  end function pio_get_ranked_filename

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    check_options
  !> @brief        check options in control file
  !! @authors      IY
  !! @param[in]    molecule     : molecule information
  !! @param[in]    boundary     : boundary information
  !! @param[in]    domain       : domain information
  !! @param[in]    trj_list     : trajectory file list information
  !! @param[in]    option       : option information
  !! @param[in]    hbond_option : hbond_option information
  !! @param[in]    totstep      : the frame number of the snapshot
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine check_options(molecule, boundary, domain, trj_list, option, &
                           hbond_option, totstep)

    ! formal arguments
    type(s_molecule),        intent(in) :: molecule
    type(s_boundary),        intent(in) :: boundary
    type(s_domain),          intent(in) :: domain
    type(s_trj_list),        intent(in) :: trj_list
    type(s_option),          intent(in) :: option
    type(s_hbond_option),    intent(in) :: hbond_option
    integer,                 intent(in) :: totstep

    ! local variables
    real(wp)                 :: hb_distance, buffer
    logical                  :: find_error


    buffer      = option%buffer
    hb_distance = hbond_option%hb_distance

    find_error  = .false.

    ! check options before the analysis loop
    !
    if (totstep == 0) then
      if (option%buffer < hb_distance) then
        write(MsgOut,'(a)') &
          'Check_Options> buffer should be larger than hb_distance'
        find_error  = .true.
      end if

      if (boundary%type == BoundaryTypeNOBC) then
        if (option%determine_box == DetermineBoxTrajectory) then
          write(MsgOut,'(a)') &
            'Check_Options> "TRAJECTORY" option is available only in PBC'
          find_error  = .true.
        end if

      else if (boundary%type == BoundaryTypePBC) then
        if (option%determine_box == DetermineboxMax) then
          write(MsgOut,*) &
            'Check_Options> "MAX" option is available only in NOBC'
          find_error  = .true.
        end if
      end if

      if (option%wrap) then
        if (trj_list%trj_type /= 2) then
          write(MsgOut,'(a)') &
          'Check_Options> when you wrap molecules, trj_type must be "COOR+BOX".'
          find_error  = .true.
        end if
      end if

      if (boundary%num_domain(1) /= 0 .and. &
          boundary%num_domain(2) /= 0 .and. &
          boundary%num_domain(3) /= 0) then

        if (product(boundary%num_domain) /= nproc_country) then
          write(MsgOut,'(a)') &
            'Check_Options> the number of process must be domain_x * domain_y * domain_z '
          find_error  = .true.
        end if
      end if
    end if

    ! check buffer distance
    !
    if (buffer < hb_distance) then
      write(MsgOut,'(a)') &
        'Check_Options> buffer distance should be larger than hb_distance.'
      find_error  = .true.
    end if

    if (find_error) then
      write(MsgOut,'(a, I8)')'Check_Options> error detected at step ',totstep
      call error_msg(&
        'Check_Options> analysis terminated because of option error')
    end if

    return

  end subroutine check_options

end module hbond_analyze_mod
