!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   qg_extract_mod
!> @brief   extract frames from trajectory files
!! @authors Norio Takase (NT), Yuji Sugita (YS), Kenta YAMADA (KYMD),
!!          Kiyoshi Yagi (KY)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module qg_extract_mod

  use qg_option_mod
  use qg_option_str_mod
  use fitting_mod
  use fitting_str_mod
  use trajectory_str_mod
  use output_str_mod
  use select_atoms_mod
  use molecules_mod
  use molecules_str_mod
  use constants_mod
  use measure_mod
  use fileio_trj_mod
  use fileio_mod
  use messages_mod
  use string_mod
  use fileio_crd_mod
  use fileio_psf_mod
  use fileio_pdb_mod
  use fileio_rst_mod
 
  implicit none
  private

  ! subroutines
  public  :: extframe
  private :: wrap_qmmm
  private :: output_snapshot
  private :: filename

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    extframe
  !> @brief        extract frame from trajectory files
  !! @authors      NT, KYMD
  !! @param[inout] molecule   : molecule information
  !! @param[inout] trj_list   : trajectory list information
  !! @param[inout] trajectory : trajectory information
  !! @param[inout] fitting    : fitting information
  !! @param[inout] option     : option information
  !! @param[inout] output     : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine extframe(molecule,   &
                      ref,        &
                      trj_list,   &
                      trajectory, &
                      fitting,    &
                      option,     &
                      output)

    ! formal arguments
    type(s_molecule),        intent(inout) :: molecule
    type(s_pdb),             intent(inout) :: ref
    type(s_trj_list),        intent(inout) :: trj_list
    type(s_trajectory),      intent(inout) :: trajectory
    type(s_fitting),         intent(inout) :: fitting
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_trj_file)         :: trj_in, trj_out
    integer                  :: iframe, irun, itrj
    integer                  :: rms_out, trr_out
    integer                  :: natom_org, alloc_stat, dealloc_stat
    real(wp), allocatable    :: ref_coord(:,:)
    type(s_rst)              :: rst
    type(s_crd)              :: crd
    type(s_pdb)              :: pdb
    logical                  :: flag_rst = .false., flag_crd = .false.
    logical                  :: flag_pdb = .false., flag_dcd = .false.
    logical                  :: first_prtqm = .true.


    ! check-only
    if (option%check_only) &
      return

    ! allocate ref_coord to save molecule%atom_coord
    !
    natom_org = molecule%num_atoms
    allocate(ref_coord(1: 3, natom_org), stat = alloc_stat)
    if (alloc_stat /= 0) call error_msg_alloc

    ! conserve molecule%atom_coord for fitting
    !
    ref_coord(1: 3, 1: natom_org) = molecule%atom_coord(1: 3, 1: natom_org)

    !
    ! for printing QM-region structure
    first_prtqm = len_trim(option%qm_atom_exp) /= 0

    iframe = 0
    do irun = 1, size(trj_list%md_steps)

      if(get_extension(trj_list%filenames(irun)) == 'dcd') then
        call open_trj(trj_in,                   &
                      trj_list%filenames(irun), &
                      trj_list%trj_format,      &
                      trj_list%trj_type,        &
                      IOFileInput)

        flag_dcd = .true.
        flag_rst = .false.
        flag_crd = .false.
        flag_pdb = .false.

      else
        trj_list%md_steps(irun)    = 1
        trj_list%ana_periods(irun) = 1
        if(get_extension(trj_list%filenames(irun)) == 'rst' .or. &
           get_extension(trj_list%filenames(irun)) == 'rsa') then
          flag_dcd = .false.
          flag_rst = .true.
          flag_crd = .false.
          flag_pdb = .false.

        else if(get_extension(trj_list%filenames(irun)) == 'crd' .or. &
                get_extension(trj_list%filenames(irun)) == 'cor') then
          flag_dcd = .false.
          flag_rst = .false.
          flag_crd = .true.
          flag_pdb = .false.

        else if(get_extension(trj_list%filenames(irun)) == 'pdb') then
          flag_dcd = .false.
          flag_rst = .false.
          flag_crd = .false.
          flag_pdb = .true.

        else
          call error_msg('extframe> &
                         &Cannot read trjfile with an unsupported file type')

        end if
      end if

      do itrj = 1, trj_list%md_steps(irun)

        ! input trj
        !
        if (flag_dcd) then
          call read_trj(trj_in, trajectory)

        else if (flag_rst) then
          call input_rst(trj_list%filenames(irun), rst)
           trajectory%pbc_box(1,1) = rst%box_size_x
           trajectory%pbc_box(2,2) = rst%box_size_y
           trajectory%pbc_box(3,3) = rst%box_size_z
           if(natom_org /= rst%num_atoms) call error_msg('extframe> &
              &Num_atoms are different between input and trjfile(RST)')
           trajectory%coord(1: 3, 1: natom_org) = rst%coord(1: 3, 1: natom_org)

        else if (flag_crd) then
          call input_crd(trj_list%filenames(irun), crd)
          trajectory%pbc_box(1,1) = 0.0_wp
          trajectory%pbc_box(2,2) = 0.0_wp
          trajectory%pbc_box(3,3) = 0.0_wp
          if(natom_org /= crd%num_atoms) call error_msg('extframe> &
                       &Num_atoms are different between input and trjfile(CRD)')

          trajectory%coord(1: 3, 1: natom_org) &
                                            = crd%atom_coord(1: 3, 1: natom_org)

        else if (flag_pdb) then
          call input_pdb(trj_list%filenames(irun), pdb)
          trajectory%pbc_box(1,1) = pdb%pbc_box(1, 1)
          trajectory%pbc_box(2,2) = pdb%pbc_box(2, 2)
          trajectory%pbc_box(3,3) = pdb%pbc_box(3, 3)
          if(natom_org /= pdb%num_atoms) call error_msg('extframe> &
                       &Num_atoms are different between input and trjfile(PDB)')

          trajectory%coord(1: 3, 1: natom_org) &
                                            = pdb%atom_coord(1: 3, 1: natom_org)

        end if

        if (mod(itrj, trj_list%ana_periods(irun)) == 0) then

          iframe = iframe + 1
          if(option%frame_no(iframe) /= 1 .and. flag_dcd) cycle 

          write(MsgOut,'(50("-"))')
          write(MsgOut,'("Frame number = ",i0,/)') iframe          


          ! wrap MM subsystems
          !
          call wrap_qmmm(molecule,   &
                         trajectory, &
                         option)


          ! substitute coordinate for select_atom and output_snapshot
          !
          molecule%atom_coord(1: 3, 1: natom_org)                              &
                                          = trajectory%coord(1: 3, 1: natom_org)


          ! selection
          !
          if(option%reconcile_num_atoms) then
            call reselect_atom(molecule,             &
                               option%qmmm_atom_exp, &
                               trajectory%coord,     &
                               option%qmmm_atom,     &
                               option%qmmm_atom_trj, &
                               tool_name='qg')
          else
            call select_atom(molecule,             &
                             option%qmmm_atom_exp, &
                             option%qmmm_atom_trj)
            write(MsgOut,'(a,i7)') 'Select_Atom> number of new selection : ', &
                                   size(option%qmmm_atom_trj%idx)
            write(MsgOut,*)

          end if

          !! wrap MM subsystems
          !!
          !call wrap_qmmm(molecule,   &
          !               trajectory, &
          !               option)

          ! fitting
          !
          if (.not.fitting%mass_weight) then
            call run_fitting(fitting,          &
                             ref_coord,        &
                             trajectory%coord, &
                             molecule%atom_coord)

          else
            call run_fitting(fitting,             &
                             ref_coord,           &
                             trajectory%coord,    &
                             molecule%atom_coord, &
                             molecule%mass)

          end if

          call output_snapshot(iframe, first_prtqm, molecule, ref, &
                               option, output)

        end if

      end do

      if (flag_dcd) then
        call close_trj(trj_in)

      else if(flag_rst) then
        call dealloc_rst_all(rst)

      else if(flag_crd) then
        call dealloc_crd_all(crd)

      else if(flag_pdb) then
        call dealloc_pdb_all(pdb)

      end if

    end do

    deallocate(ref_coord, stat = dealloc_stat)
    if (dealloc_stat /= 0) call error_msg_dealloc

    return

  end subroutine extframe

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    wrap_qmmm
  !> @brief        wrap solvent with the COM of QM subsystem centered (nonBC)
  !! @authors      KYMD
  !! @param[in]    molecule   : molecule information
  !! @param[inout] trajectory : trajectory information
  !! @param[in]    option     : option information
  !! @note         setup_pbc_correct must be performed previously
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine wrap_qmmm(molecule, trajectory, option)

    ! formal argments
    type(s_molecule),           intent(in)    :: molecule
    type(s_trajectory), target, intent(inout) :: trajectory
    type(s_option),             intent(in)    :: option

    ! local variables
    real(wp)                    :: origin(3), com(3)
    real(wp)                    :: boxx, boxy, boxz
    real(wp)                    :: boxx_inv, boxy_inv, boxz_inv
    integer                     :: natom, nmole
    integer                     :: i, j, initial, final

    real(wp),           pointer :: coord(:,:)


    coord => trajectory%coord

    natom  = molecule%num_atoms
    nmole  = molecule%num_molecules
    origin = compute_com(coord, molecule%mass, option%origin_atom%idx)

    boxx     = trajectory%pbc_box(1,1)
    boxy     = trajectory%pbc_box(2,2)
    boxz     = trajectory%pbc_box(3,3)
    if(boxx * boxy * boxz < EPS) then
      write(MsgOut, '(9x, a)') '* Skip wrapping solvent molecules'
      return
    end if
    boxx_inv = 1.0_wp / boxx
    boxy_inv = 1.0_wp / boxy
    boxz_inv = 1.0_wp / boxz

    do i = 1, nmole

      initial  = molecule%molecule_atom_no(i)
      if (i /= nmole) then
        final = molecule%molecule_atom_no(i + 1) - 1
      else
        final = natom
      end if

      ! compute center of mass of a molecule
      !
      com(1) = 0.0_wp
      com(2) = 0.0_wp
      com(3) = 0.0_wp

      do j = initial, final
        coord(1, j) = coord(1, j) - origin(1)
        com(1)      = com(1) + coord(1, j) * molecule%mass(j)
        coord(2, j) = coord(2, j) - origin(2)
        com(2)      = com(2) + coord(2, j) * molecule%mass(j)
        coord(3, j) = coord(3, j) - origin(3)
        com(3)      = com(3) + coord(3, j) * molecule%mass(j)
      end do

      com(1) = com(1) / molecule%molecule_mass(i)
      com(2) = com(2) / molecule%molecule_mass(i)
      com(3) = com(3) / molecule%molecule_mass(i)

      ! move molecule into the unit cell
      !
      do j = initial, final
        coord(1, j) = coord(1, j) - boxx * anint(com(1)*boxx_inv) + origin(1)
        coord(2, j) = coord(2, j) - boxy * anint(com(2)*boxy_inv) + origin(2)
        coord(3, j) = coord(3, j) - boxz * anint(com(3)*boxz_inv) + origin(3)
      end do

    end do

    return

  end subroutine wrap_qmmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_snapshot
  !> @brief        save frame to snapshot file
  !! @authors      KYMD, KY
  !! @param[in   ] iframe      : structure index
  !! @param[inout] first_prtqm : flag of printing QM region
  !! @param[inout] molecule    : molecule information
  !! @param[inout] ref        : pdb information
  !! @param[inout] option      : option information
  !! @param[inout] output      : output information
  !! @note         this subroutine was originally made by MT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine output_snapshot(iframe, first_prtqm, molecule, ref, option, output)

    ! formal arguments
    integer,                 intent(in   ) :: iframe
    logical,                 intent(inout) :: first_prtqm
    type(s_molecule),        intent(inout) :: molecule
    type(s_pdb),             intent(inout) :: ref 
    type(s_option),          intent(inout) :: option
    type(s_output),          intent(inout) :: output

    ! local variables
    type(s_crd)              :: tmp_crd
    type(s_psf)              :: tmp_psf
    type(s_pdb)              :: tmp_pdb
    character(MaxFilename)   :: line
    integer                  :: ipos_dot
    logical                  :: ioext


    line = ''
    ioext = .true.

    ! Punch out COORDINATE file
    !
    if (len_trim(output%qmmm_crdfile) /= 0) then

      select case (option%coord_format)

      case (CrdFormatCharmm)

        call export_molecules(molecule, option%qmmm_atom_trj, crd=tmp_crd) 

        call output_crd(filename(output%qmmm_crdfile, iframe), tmp_crd, ioext)
        call dealloc_crd_all(tmp_crd)

        write(line, '(a)') trim(filename(output%qmmm_crdfile, iframe))

      case default

        call error_msg('output_snapshot> &
                       &Formats other than CHARMM_CARD not implemented yet')

      end select

    endif

    ! Punch out PDB file for analysis 
    !
    if (len_trim(output%qmmm_pdbfile) /= 0) then

      call export_molecules(molecule, option%qmmm_atom_trj, pdb=tmp_pdb)
      tmp_pdb%ter_rec   = ref%ter_rec
      call output_pdb(filename(output%qmmm_pdbfile, iframe), tmp_pdb)
      call dealloc_pdb_all(tmp_pdb)

      if (len_trim(line) == 0)  &
         write(line, '(a)') trim(filename(output%qmmm_pdbfile, iframe))

    end if

    ! Punch out PSF 
    !
    if (len_trim(output%qmmm_psffile) /= 0) then

      ! Duplicate psf
      call duplicate_psf(option%dup_psf, tmp_psf)

      call export_molecules(molecule, option%qmmm_atom_trj, psf=tmp_psf)
      call output_psf(filename(output%qmmm_psffile, iframe), tmp_psf)
      call dealloc_psf_all(tmp_psf)

      if (len_trim(line) == 0)  &
         write(line, '(a)') trim(filename(output%qmmm_psffile, iframe))

    end if

    write(MsgOut, '(9x, 2a)') '* Save this frame to ', trim(line)
    write(MsgOut,*)

    ! Print QM-region structure
    !
    if(first_prtqm) then

      ipos_dot = index(line, '.', back = .true.) - 1
      if (ipos_dot < 0) ipos_dot = len_trim(line)
      if (ipos_dot + 13 < MaxFilename) then
        write(line, '(a)') line(1: ipos_dot)//"_qmregion.pdb"

      else
        line = ''
        write(line, '(a)') "qmregion.pdb"

      endif

      call export_molecules(molecule, option%qm_atom, pdb=tmp_pdb) 
      tmp_pdb%ter_rec   = ref%ter_rec
      call output_pdb(line, tmp_pdb)
      call dealloc_pdb_all(tmp_pdb)

      write(MsgOut, '(9x, 2a)') '* Save QM-region structure in this frame to ', trim(line)
      write(MsgOut,*)

      first_prtqm = .false.

    endif

    return

  end subroutine output_snapshot

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    filename
  !> @brief        insert frame index into {} in the basename
  !! @authors      KYMD
  !! @param[in]    basename      : basename of filename
  !! @param[in]    num           : index
  !! @note         this subroutine was originally made by NT and TM
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function filename(basename, num)

    ! return
    character(MaxFilename)   :: filename

    ! formal arguments
    character(*),            intent(in)    :: basename
    integer,                 intent(in)    :: num

    ! local variables
    integer                  :: bl, br
    !character(100)           :: fid


    bl = index(basename, '{', back=.true.)
    br = index(basename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) then
      write(MsgOut, '(a)') ' (Warning) {} is not found in  '//trim(basename)
      filename=basename
    else
      !write(MsgOut, '(2a, i0, a)') '         * Save this frame to ', &
      !  basename(:bl-1), num, trim(basename(br+1:))
      write(filename, '(a, i0, a)') basename(:bl-1), num, trim(basename(br+1:))
    end if

    return

  end function filename

end module qg_extract_mod
