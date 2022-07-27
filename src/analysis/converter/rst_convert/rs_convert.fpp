!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   rs_convert_mod
!> @brief   convert restart files
!! @authors Norio Takase (NT), Chigusa Kobayashi (CK)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../../config.h"
#endif

module rs_convert_mod

  use rs_option_mod
  use rs_option_str_mod
  use output_str_mod
  use input_str_mod
  use fileio_mod
  use fileio_rst_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_ambcrd_mod
  use fileio_namd_xyz_mod
  use fileio_namd_xsc_mod
  use messages_mod
  use molecules_str_mod
  use constants_mod
  use string_mod
  use pbc_correct_mod
  use trajectory_str_mod
 
  implicit none
  private

  ! parameters
  character(*), private, parameter :: RstFileTypes(3) = (/'MIN ','MD  ','REMD'/)

  ! subroutines
  public  :: convert
  private :: remd_rst_convert
  private :: rst2various
  private :: input_namd
  private :: convert_pdb
  private :: convert_charmm
  private :: convert_namd
  private :: convert_ambcrd
  private :: get_filename

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert
  !> @brief        convert restart file
  !! @authors      NT
  !! @param[in]    input    : input information
  !! @param[in]    molecule : molecule information
  !! @param[in]    option   : option information
  !! @param[in]    output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine convert(input, molecule, option, output)

    ! formal arguments
    type(s_input),    intent(in) :: input
    type(s_molecule), intent(in) :: molecule
    type(s_option),   intent(in) :: option
    type(s_output),   intent(in) :: output

    ! local variables
    type(s_rst)                  :: rst
    integer                      :: bl, br

    ! check check-only
    if (option%check_only) &
      return

    ! open Source restart file and convert
    select case (option%rstin_format)

    case (RstFormatNAMD)
      call input_namd(input, option, rst)

    case (RstFormatGENESIS)

      bl = index(input%rstfile, '{', back=.true.)
      br = index(input%rstfile, '}', back=.true.)


      ! rstout_type = 'AUTO'
      if (option%rstout_type == RstTypeAuto) then
        if ((bl /= 0) .and. (br /= 0) .and. (bl < br)) then
          call remd_rst_convert(input, molecule, option, output)
          return

        else
          call input_rst(input%rstfile, rst)

        end if

      ! rstout_type = 'REMD'
      else if (option%rstout_type == RstTypeREMD) then
        if ((bl /= 0) .and. (br /= 0) .and. (bl < br)) then
          call remd_rst_convert(input, molecule, option, output)
          return

        else
          call error_msg('When you want to convert REMD rst file, you need to use {}')
        end if

      ! rstout_type = 'MD' or 'MIN'
      else
        call input_rst(input%rstfile, rst)

      end if

    end select

    call rst2various(rst, molecule, option, output)

    ! deallocate memory
    call dealloc_rst_all(rst)

    write(MsgOut,'(A)') ' '

    return

  end subroutine convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    remd_rst_convert
  !> @brief        convert restart file (For REMD)
  !! @authors      Daisuke Matsuoka (DM)
  !! @param[in]    input    : input information
  !! @param[in]    molecule : molecule information
  !! @param[in]    option   : option information
  !! @param[in]    output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine remd_rst_convert(input, molecule, option, output)

    ! formal arguments
    type(s_input),    intent(in) :: input
    type(s_molecule), intent(in) :: molecule
    type(s_option),   intent(in) :: option
    type(s_output),   intent(in) :: output

    ! local variables
    type(s_rst)                  :: rst
    type(s_output)               :: output_replica
    integer                      :: nreplicas, i_repid
    character(MaxFileName)       :: infil, outfil

    ! setup
    infil = get_filename(input%rstfile, 1)
    call input_rst(infil, rst)
    nreplicas = rst%nreplicas

    if (option%nreplicas > 0) then
      nreplicas = option%nreplicas
    end if

    ! convert rst file for each replica
    do i_repid = 1, nreplicas
      infil  = get_filename(input%rstfile, i_repid)

      write(MsgOut,'(A)') ' '
      write(MsgOut,'(A)') 'Read_Rst_FileName> '//trim(infil)
      call input_rst(infil, rst)

      select case (option%rstout_format)

      case (RstFormatGENESIS)
        output_replica%rstfile = get_filename(output%rstfile, i_repid)

      case (RstFormatPDB)
        output_replica%pdbfile = get_filename(output%pdbfile, i_repid)

      case (RstFormatCHARMM)
        output_replica%crdfile = get_filename(output%crdfile, i_repid-1)

      case (RstFormatAMBER)
        output_replica%ambcrdfile = get_filename(output%ambcrdfile, i_repid)

      end select

      call rst2various(rst, molecule, option, output_replica)

      ! deallocate memory
      call dealloc_rst_all(rst)
    end do

    write(MsgOut,'(A)') ' '

    return

  end subroutine remd_rst_convert

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    rst2various
  !> @brief        convert rst format to various formats
  !! @authors      DM
  !! @param[in]    rst      : GENESIS restart information
  !! @param[in]    molecule : molecule information
  !! @param[in]    option   : option information
  !! @param[in]    outfil   : output file name
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine rst2various(rst, molecule, option, output)

    ! formal arguments
    type(s_rst),      intent(inout) :: rst
    type(s_molecule), intent(in)    :: molecule
    type(s_option),   intent(in)    :: option
    type(s_output),   intent(in)    :: output

    ! local variables
    type(s_trajectory)       :: trajectory
    integer                  :: iatm

    ! wrap_molecule
    call alloc_trajectory(trajectory, rst%num_atoms)
    do iatm = 1, size(molecule%atom_no)
      trajectory%coord(1:3, iatm) = rst%coord(1:3, iatm) + option%shift_coord(1:3)
    end do

    trajectory%pbc_box(1,1) = rst%box_size_x
    trajectory%pbc_box(2,2) = rst%box_size_y
    trajectory%pbc_box(3,3) = rst%box_size_z

    call run_pbc_correct(option%pbcc_mode, molecule, trajectory)

    rst%coord = trajectory%coord
    call dealloc_trajectory(trajectory)

    ! write restart file
    select case (option%rstout_format)

    case (RstFormatNAMD)
      call convert_namd(output%coorfile, output%velfile,  output%xscfile, &
                       option, rst)

    case (RstFormatGENESIS)
      call output_rst(output%rstfile, rst)

    case (RstFormatPDB)
      call convert_pdb(output%pdbfile, molecule, rst)

    case (RstFormatCHARMM)
      call convert_charmm(output%crdfile, molecule, rst)

    case (RstFormatAMBER)
      call convert_ambcrd(output%ambcrdfile, rst)

    end select

    write(MsgOut,'(A)') ' '

    return

  end subroutine rst2various

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_namd
  !> @brief        read NAMD restart file
  !! @authors      NT
  !! @param[in]    input  : input information
  !! @param[in]    option : option information
  !! @param[inout] rst    : GENESIS restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_namd(input, option, rst)

    ! formal arguments
    type(s_input),           intent(in)    :: input
    type(s_option),          intent(in)    :: option
    type(s_rst),             intent(inout) :: rst

    !local variables
    type(s_namd_xyz)         :: crd, vel
    type(s_namd_xsc)         :: xsc


    ! input NAMD restart informations
    call input_namd_xyz(input%coorfile, crd)
    call input_namd_xsc(input%xscfile, xsc)

    if (len_trim(input%velfile) > 0) &
      call input_namd_xyz(input%velfile, vel)


    ! rstfile type
    select case(option%rstout_type)

    case (RstTypeAuto)

      if (.not. allocated(vel%v)) then
        rst%rstfile_type = RstfileTypeMin
      else
        rst%rstfile_type = RstfileTypeMd
      end if

    case (RstTypeMin)

      rst%rstfile_type = RstfileTypeMin

    case (RstTypeMD)

      rst%rstfile_type = RstfileTypeMd

    end select

    write(MsgOut,'(A)') 'Run_Namd_Rstcnv> Genesis restart file type: '// &
         RstFileTypes(rst%rstfile_type)

    ! check coordinates and velocities
    if (rst%rstfile_type == rstfileTypeMd) then

      if (.not. allocated(vel%v)) &
        call error_msg('Run_Namd_Rstcnv> velocity file is not indicated.')

      if (size(crd%v) /= size(vel%v)) &
        call error_msg('Run_Namd_Rstcnv> velocity atom count is different'//&
        ' from coordinates.')

    end if

    rst%num_atoms = size(crd%v(1,:))
    call alloc_rst(rst, RestartAtom, rst%num_atoms)


    ! box size
    rst%box_size_x = xsc%a(1)
    rst%box_size_y = xsc%b(2)
    rst%box_size_z = xsc%c(3)

    write(MsgOut,'(A,3F8.3)')  'Run_Namd_Rstcnv> box size : ', &
                               rst%box_size_x, &
                               rst%box_size_y, &
                               rst%box_size_z

    ! min: energy, delta_r
    rst%energy  = option%min_energy
    rst%delta_r = option%min_delta_r

    if (rst%rstfile_type == rstfileTypeMin) then
      write(MsgOut,'(A,F8.3)') 'Run_Namd_Rstcnv> energy   : ', rst%energy
      write(MsgOut,'(A,F8.3)') 'Run_Namd_Rstcnv> delta_r  : ', rst%delta_r
    end if

    ! md : iseed, ndof, thermo_moment, baro_moment
    rst%iseed               = option%md_iseed
    rst%thermostat_momentum = option%md_thermo_moment
    rst%barostat_momentum   = option%md_baro_moment

    if (rst%rstfile_type == rstfileTypeMd) then
      write(MsgOut,'(A,I16)')   'Run_Namd_Rstcnv> iseed    : ', rst%iseed
      write(MsgOut,'(A,F16.9)') 'Run_Namd_Rstcnv> thermostat-momentum : ', &
                               rst%thermostat_momentum
      write(MsgOut,'(A,F16.9)')'Run_Namd_Rstcnv> barostat-momentum   : ', &
                               rst%barostat_momentum(1)
      write(MsgOut,'(A,F16.9)')'                                       ', &
                               rst%barostat_momentum(2)
      write(MsgOut,'(A,F16.9)')'                                       ', &
                               rst%barostat_momentum(3)
    end if

    ! coordinate
    rst%coord(:,:) = crd%v(:,:)

    ! velocity
    if (rst%rstfile_type == rstfileTypeMd) &
      rst%velocity(:,:) = vel%v(:,:)

    write(MsgOut,'(A)') ' '


    ! deallocate memory
    call dealloc_namd_xyz(vel)
    call dealloc_namd_xyz(crd)

    return

900 call error_msg('Run_Rstcnv> Input file open error.')

  end subroutine input_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_pdb
  !> @brief        convert pdb coordinate file
  !! @authors      DM
  !! @param[in]    pdbfile  : PDB file name
  !! @param[in]    molecule : molecule information
  !! @param[inout] rst      : GENESIS restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_pdb(pdbfile, molecule, rst)

    ! formal arguments
    character(*),            intent(in)    :: pdbfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_rst),             intent(inout) :: rst

    ! local variable
    type(s_pdb)                            :: pdb

    call init_pdb(pdb)
    pdb%num_atoms = rst%num_atoms

    if (pdb%num_atoms /= size(molecule%atom_no)) &
      call error_msg('OUTPUT_CHARMM> number of atoms in rst is not equal to '// &
                                    'number of atoms in psf/prmtop/grotop file')

    call alloc_pdb(pdb, PdbAtom, pdb%num_atoms)

    pdb%cryst_rec    = .true.
    pdb%hetatm_rec   = .true.
    pdb%end_rec      = .true.
    pdb%segment      = .true.

    pdb%atom_no      = molecule%atom_no
    pdb%atom_name    = molecule%atom_name
    pdb%residue_no   = molecule%residue_no
    pdb%residue_name = molecule%residue_name
    pdb%chain_id     = molecule%chain_id
    pdb%segment_name = molecule%segment_name
    pdb%atom_coord   = rst%coord
    pdb%pbc_box(1,1) = rst%box_size_x
    pdb%pbc_box(2,2) = rst%box_size_y
    pdb%pbc_box(3,3) = rst%box_size_z

    call output_pdb(pdbfile, pdb)

    call dealloc_pdb_all(pdb)

    return

  end subroutine convert_pdb

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_charmm
  !> @brief        convert CHARMM coordinate file
  !! @authors      DM
  !! @param[in]    crdfile  : CHARMM crd file name
  !! @param[in]    molecule : molecule information
  !! @param[inout] rst      : GENESIS restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine convert_charmm(crdfile, molecule, rst)

    ! formal arguments
    character(*),            intent(in)    :: crdfile
    type(s_molecule),        intent(in)    :: molecule
    type(s_rst),             intent(inout) :: rst

    ! local variable
    type(s_crd)                            :: charmm_crd

    call init_crd(charmm_crd)
    charmm_crd%num_atoms = rst%num_atoms

    if (charmm_crd%num_atoms /= size(molecule%atom_no)) &
      call error_msg('OUTPUT_CHARMM> number of atoms in rst is not equal to '// &
                                    'number of atoms in psf/prmtop/grotop file')

    call alloc_crd(charmm_crd, CrdAtom, charmm_crd%num_atoms)

    charmm_crd%atom_no      = molecule%atom_no
    charmm_crd%atom_name    = molecule%atom_name
    charmm_crd%residue_no   = molecule%residue_no
    charmm_crd%residue_name = molecule%residue_name
    charmm_crd%segment_name = molecule%segment_name
    charmm_crd%atom_coord   = rst%coord

    write(MsgOut,'(A,3F8.3)')  'Run_CHARMM_convert> box size : ', &
                               rst%box_size_x, &
                               rst%box_size_y, &
                               rst%box_size_z

    call output_crd(crdfile, charmm_crd, .true.)

    call dealloc_crd_all(charmm_crd)

    return

  end subroutine convert_charmm

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    output_namd
  !> @brief        a driver subroutine for writing NAMD restart file
  !! @authors      CK
  !! @param[in]    coorfile : filename of coorfile
  !! @param[in]    velfile  : filename of velocity file
  !! @param[in]    xscfile  : extended system information file
  !! @param[in]    option   : option information
  !! @param[inout] rst      : GENESIS restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_namd(coorfile, velfile, xscfile, option, rst)

    ! formal arguments
    character(*),   intent(in)    :: coorfile
    character(*),   intent(in)    :: velfile
    character(*),   intent(in)    :: xscfile
    type(s_option), intent(in)    :: option
    type(s_rst),    intent(inout) :: rst

    ! local variables
    type(s_namd_xyz)            :: crd, vel
    type(s_namd_xsc)            :: xsc

    if (rst%rstfile_type == RstfileTypeRemd)  &
       call error_msg('Convert_Genesis> Remd type restart is not allowed')

    ! initialize coor, vel, xsc
    call alloc_namd_xyz(crd, rst%num_atoms)
    call alloc_namd_xyz(vel, rst%num_atoms)

    xsc%step = option%step
    xsc%a(:) = 0.0_dp
    xsc%b(:) = 0.0_dp
    xsc%c(:) = 0.0_dp
    xsc%o(:) = 0.0_dp
    xsc%s(:) = 0.0_dp

    crd%v(:,:) = 0.0_dp
    vel%v(:,:) = 0.0_dp

    ! box size
    xsc%a(1) = rst%box_size_x
    xsc%b(2) = rst%box_size_y
    xsc%c(3) = rst%box_size_z
    xsc%o(1) = option%origin_x
    xsc%o(2) = option%origin_y
    xsc%o(3) = option%origin_z

    write(MsgOut,'(A,3F8.3)')  'Run_Namd_Rstcnv> box size : ', &
                               xsc%a(1), &
                               xsc%b(2), &
                               xsc%c(3)
    write(MsgOut,'(A,3F8.3)')  'Run_Namd_Rstcnv> cell origin : ', &
                               xsc%o(1:3)

    ! coordinate
    crd%v(:,:) = dble(rst%coord(:,:))

    ! velocity
    vel%v(:,:) = dble(rst%velocity(:,:))

    ! open xscfile
    !
    call output_namd_xsc(xscfile, xsc)

    ! open coorfile
    !
    call output_namd_xyz(coorfile, crd)

    ! open velfile
    !
    call output_namd_xyz(velfile, vel)

    ! deallocate memory
    call dealloc_namd_xyz(crd)
    call dealloc_namd_xyz(vel)

     return

  end subroutine convert_namd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    convert_ambcrd
  !> @brief        convert AMBER coordinate file
  !! @authors      DM
  !! @param[in]    abmcrdfile : ABMER crd file name
  !! @param[inout] rst        : GENESIS restart information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine convert_ambcrd(ambcrdfile, rst)

    ! formal arguments
    character(*),            intent(in)    :: ambcrdfile
    type(s_rst),             intent(inout) :: rst

    ! local variable
    type(s_ambcrd)                         :: ambcrd


    call init_ambcrd(ambcrd)
    ambcrd%title = 'CONVERTED BY RST_CONVERT IN GENESIS ANALYSIS TOOLS'

    ambcrd%num_atoms    = rst%num_atoms
    ambcrd%box_rec      = .true.
    ambcrd%box_angl_rec = .true.
    ambcrd%box_size(1)  = rst%box_size_x
    ambcrd%box_size(2)  = rst%box_size_y
    ambcrd%box_size(3)  = rst%box_size_z
    ambcrd%box_angl     = 90.0_wp

    if ((rst%rstfile_type == RstfileTypeMD) .or. &
        (rst%rstfile_type == RstfileTypeREMD)) then

      ambcrd%velocity_rec = .true.
    end if

    call alloc_ambcrd(ambcrd, PrmcrdAtom, ambcrd%num_atoms)

    ambcrd%atom_coord    = rst%coord
    ambcrd%atom_velocity = rst%velocity

    call output_ambcrd(ambcrdfile, ambcrd)

    call dealloc_ambcrd_all(ambcrd)

    return

  end subroutine convert_ambcrd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    get_filename
  !> @brief        insert snapshot index into {} in the filename
  !! @authors      TM
  !! @param[in]    filename      : filename
  !! @param[in]    no            : index
  !! @note         this subroutine was originally made by NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  function get_filename(filename, no)

    ! return
    character(Maxfilename)   :: get_filename

    ! formal arguments
    character(*),            intent(in)    :: filename
    integer,                 intent(in)    :: no

    ! local variables
    integer                  :: bl, br
    character(100)           :: fid


    bl = index(filename, '{', back=.true.)
    br = index(filename, '}', back=.true.)

    if (bl == 0 .or. br == 0 .or. bl > br) &
      call error_msg('Get_Filename> {} is not found in the output trjfile name')

    write(fid,'(i0)') no
    get_filename = filename(:bl-1) // trim(fid) // filename(br+1:)

    return

  end function get_filename

end module rs_convert_mod
