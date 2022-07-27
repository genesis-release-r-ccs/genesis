!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_control_mod
!> @brief   read parameters and data for Parallel I/O restart setup
!! @authors Norio Takase (NT), Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module pr_control_mod

  use sp_dynamics_mod
  use sp_output_mod
  use sp_input_mod
  use sp_ensemble_mod
  use sp_restraints_mod
  use sp_constraints_mod
  use sp_boundary_mod
  use sp_energy_mod
  use sp_ensemble_str_mod
  use select_mod
  use fileio_control_mod
  use messages_mod
  use mpi_parallel_mod
  use constants_mod

  implicit none
  private

  ! structures
  type, public :: s_optional_info

    ! data for section output
    character(200)         :: cachepath   = './cache'
    logical                :: keep_cache  = .false.
    logical                :: convert     = .false.
    logical                :: set_restart = .false.

    ! data for section boundary
    integer                :: domain_xyz = 0

  end type s_optional_info

  type, public :: s_ctrl_data

    ! data for section input
    type(s_inp_info)       :: inp_info

    ! data for section output
    type(s_out_info)       :: out_info

    ! data for section energy
    type(s_ene_info)       :: ene_info

    ! data for section dynamics
    type(s_dyn_info)       :: dyn_info

    ! data for section constraints
    type(s_cons_info)      :: cons_info

    ! data for section ensemble
    type(s_ens_info)       :: ens_info

    ! data for section boundary
    type(s_boundary_info)  :: bou_info

    ! data for section selection
    type(s_sel_info)       :: sel_info

    ! data for section restraint
    type(s_res_info)       :: res_info

    ! additional information
    type(s_optional_info)  :: opt_info

  end type s_ctrl_data

  ! subroutines
  public  :: usage
  public  :: control
  private :: init_ctrl_data
  private :: read_ctrl_output2
  private :: read_ctrl_dynamics2
  private :: read_ctrl_boundary2

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    usage
  !> @brief        show usage of prst_setup
  !! @authors      NT, TM
  !! @param[in]    arg1 : arguments in execution
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine usage(arg1)

    ! formal arguments
    character(*),            intent(out) :: arg1

    ! local variables
    character(80)            :: arg2
    integer                  :: iargc


    call getarg(1, arg1)

    if (iargc() < 1      .or. &
         arg1 == '-help' .or. &
         arg1 == '-h'    .or. &
         arg1 == '-HELP' .or. &
         arg1 == '-H') then

      if (main_rank) then

        call getarg(2, arg2)

        if (iargc() < 2 .or. arg2 /= 'ctrl') then
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# normal usage'
          write(MsgOut,'(A)') '  % ./prst_setup INPSET'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '# usage to see control parameters'
          write(MsgOut,'(A)') '  % ./prst_setup -h ctrl'
          write(MsgOut,'(A)') ' '
        else
          write(MsgOut,'(A)') '# control parameters in prst_setup'
          write(MsgOut,'(A)') ' '
          write(MsgOut,'(A)') '[INPUT]'
          write(MsgOut,'(A)') '##  CHARMM Force Field'
          write(MsgOut,'(A)') 'topfile       = sample.top      # CHARMM topology file'
          write(MsgOut,'(A)') 'parfile       = sample.par      # CHARMM parameter file'
          write(MsgOut,'(A)') 'psffile       = sample.psf      # CHARMM protein structure file'
          write(MsgOut,'(A)') 'pdbfile       = sample.pdb      # PDB coordinate file'
          write(MsgOut,'(A)') '# crdfile       = sample.crd      # CHARMM coordinates file'
          write(MsgOut,'(A)') '# reffile       = sample.pdb      # Reference PDB coordinate file'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '##  AMBER Force Field'
          write(MsgOut,'(A)') '# prmtopfile    = sample.top   # AMBER parameter topology file'
          write(MsgOut,'(A)') '# ambcrdfile    = sample.crd   # AMBER coordinate file'
          write(MsgOut,'(A)') '# ambreffile    = sample.crd   # Reference AMBER coordinate file'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '##  GROMACS Force Field'
          write(MsgOut,'(A)') '# grotopfile    = sample.grotop # GROMACS topology file'
          write(MsgOut,'(A)') '# grocrdfile    = sample.grocrd # GROMACS coordinate file'
          write(MsgOut,'(A)') '# groreffile    = sample.grocrd # Reference GROMACS coordinate file'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '# rstfile       = sample().rst    # Domain restart files'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[OUTPUT]      '
          write(MsgOut,'(A)') 'rstfile       = ./output().rst  # GENESIS restart file'
          write(MsgOut,'(A)') '# selfile       = sample.sel    # selection list file'
          write(MsgOut,'(A)') '# pdbfile       = ./output.pdb    # PDB file'
          write(MsgOut,'(A)') 'cachepath     = ./cache'
          write(MsgOut,'(A)') 'keep_cache    = NO'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[ENERGY]'
          write(MsgOut,'(A)') 'pairlistdist  = 13.5'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[DYNAMICS]'
          write(MsgOut,'(A)') '# iseed         = -1        # random number seed'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[CONSTRAINTS]'
          write(MsgOut,'(A)') 'rigid_bond    = YES         # constraints all bonds involving hydrogen'
          write(MsgOut,'(A)') '# shake_iteration = 500     # max number of SHAKE/RATTLE iterations'
          write(MsgOut,'(A)') '# shake_tolerance = 1.0e-8  # SHAKE/RATTLE tolerance (Ang)'
          write(MsgOut,'(A)') '# water_model     = TIP3    # water model'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[ENSEMBLE]'
          write(MsgOut,'(A)') 'ensemble      = NVE       # [NVE,NVT,NPT,NPAT]'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[SELECTION]'
          write(MsgOut,'(A)') '# group1        = atomname:CA  # selection group 1'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[RESTRAINTS]'
          write(MsgOut,'(A)') '# nfunctions    = 1       # number of restraint function 1'
          write(MsgOut,'(A)') '# function1     = POSI    # [DIST,POSI,RMSD] [MASS] kind of function 1'
          write(MsgOut,'(A)') '# constant1     = 10.0    # force constant(x,y,z) for function 1'
          write(MsgOut,'(A)') '# select_index1 = 1       # restraint groups in function 1'
          write(MsgOut,'(A)') ''
          write(MsgOut,'(A)') '[BOUNDARY]'
          write(MsgOut,'(A)') 'box_size_x    = 0.0       # box size (x) in [PBC]'
          write(MsgOut,'(A)') 'box_size_y    = 0.0       # box size (y) in [PBC]'
          write(MsgOut,'(A)') 'box_size_z    = 0.0       # box size (z) in [PBC]'
          write(MsgOut,'(A)') 'domain_xyz    = 8         # # of total domain'
          write(MsgOut,'(A)') '# domain_x       = 2      # # of domain for x-axis'
          write(MsgOut,'(A)') '# domain_y       = 2      # # of domain for y-axis'
          write(MsgOut,'(A)') '# domain_z       = 2      # # of domain for z-axis'
          write(MsgOut,'(A)') ' '
        end if

        stop

      end if

    else

      if (.not. main_rank) &
        return

      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*                    PRST_SETUP                    *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*      Generater of parallel I/O restart files     *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'*                Developed by RIKEN                *'
      write(MsgOut,'(A)')'*                                                  *'
      write(MsgOut,'(A)')'****************************************************'
      write(MsgOut,'(A)')' '
    end if

    return

  end subroutine usage

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    control
  !> @brief        open/read/close control files
  !! @authors      NT, TM
  !! @param[in]    filename  file name of control file
  !! @param[inout] ctrl_data information of control parameters
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine control(filename, ctrl_data)

    ! formal arguments
    character(*),            intent(in)    :: filename
    type(s_ctrl_data),       intent(inout) :: ctrl_data

    ! local variables
    integer                  :: handle


    ! open control file
    !
    call open_ctrlfile(filename, handle)

    if (handle == 0) &
      call error_msg('Control> File IO Error')

    ! initialize variables
    !
    call init_ctrl_data(ctrl_data)

    ! read input section
    !
    call read_ctrl_input (handle, ctrl_data%inp_info)

    ! read output section
    !
    call read_ctrl_output2(handle, ctrl_data%opt_info)
    call read_ctrl_output (handle, ctrl_data%out_info)

    ! read energy section
    !
    call read_ctrl_energy(handle, ctrl_data%ene_info)

    ! read dynamics section
    !
    call read_ctrl_dynamics2(handle, ctrl_data%dyn_info)

    ! read constraints section
    !
    call read_ctrl_constraints(handle, ctrl_data%cons_info)

    ! read ensemble section
    !
    call read_ctrl_ensemble(handle, ctrl_data%ens_info)

    ! read boundary section
    !
    call read_ctrl_boundary2(handle, ctrl_data%opt_info)
    call read_ctrl_boundary (handle, ctrl_data%bou_info)

    ! read selection section
    !
    call read_ctrl_selection(handle, ctrl_data%sel_info)

    ! read restraint section
    !
    call read_ctrl_restraints(handle, ctrl_data%res_info)

    ! close control file
    !
    call close_ctrlfile(handle)

    return

  end subroutine control

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    init_ctrl_data
  !> @brief        initialize ctrl_data variables
  !! @authors      NT
  !! @param[inout] ctrl_data information of control parameters [s_ctrl_data]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine init_ctrl_data(ctrl_data)

    ! formal arguments
    type(s_ctrl_data),       intent(inout) :: ctrl_data


    ! initialize energy information variables
    !

    ctrl_data%ene_info%switchdist   = 0.000000_wp
    ctrl_data%ene_info%cutoffdist   = 0.000001_wp
    ctrl_data%ene_info%dielec_const = 0.000000_wp


    return

  end subroutine init_ctrl_data

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_output2
  !> @brief        read parameters of output section
  !! @param[in]    handle   : handle for control file
  !! @param[out]   opt_info : optional information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_output2(handle, opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Output'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_optional_info),   intent(inout) :: opt_info


    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string (handle, Section, 'cachepath',  & !TODO re-name
                              opt_info%cachepath)

    call read_ctrlfile_logical(handle, Section, 'keep_cache', & !TODO delete
                              opt_info%keep_cache)

    call read_ctrlfile_logical(handle, Section, 'old_convert', &
                              opt_info%convert)

    call read_ctrlfile_logical(handle, Section, 'restart', &
                              opt_info%set_restart)

    if (main_rank) then

      write(MsgOut,'(a)')   'Read_Ctrl_Output2> '
      write(MsgOut,'(a,a)') '  Cache file path = ', trim(opt_info%cachepath)
      if (opt_info%keep_cache) then
        write(MsgOut,'(a,a)') '  Keep cache file = yes'
      else
        write(MsgOut,'(a,a)') '  Keep cache file = no'
      end if
      if (opt_info%convert) then
        write(MsgOut,'(a,a)') 'Convert old to new = yes'
      endif
      if (opt_info%set_restart) then
        write(MsgOut,'(a,a)') '  Set restart      = yes'
      endif
      write(MsgOut,'(a)')   ' '

    end if

    return

  end subroutine read_ctrl_output2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_dynamics2
  !> @brief        read DYNAMICS section in the control file
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   dyn_info : DYNAMICS section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_dynamics2(handle, dyn_info)

    ! parameters
    character(*),            parameter     :: Section = 'Dynamics'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_dyn_info),        intent(inout) :: dyn_info


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'Iseed',         &
                               dyn_info%iseed)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Dynamics> Parameters of MD simulation'

      write(MsgOut,'(A20,I10)') &
            '  iseed           = ', dyn_info%iseed
      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_dynamics2

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_boundary2
  !> @brief        read parameters of boundary section
  !! @param[in]    handle   : handle for control file
  !! @param[inout] opt_info : optional information
  !! @authors      NT
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_boundary2(handle, opt_info)

    ! parameters
    character(*),            parameter     :: Section = 'Boundary'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_optional_info),   intent(inout) :: opt_info


    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_integer(handle, Section, 'domain_xyz',  &
                               opt_info%domain_xyz)


    if (main_rank) then

      write(MsgOut,'(a)')     'Read_Ctrl_Boundary2> '
      write(MsgOut,'(a,i10)') '  Domain xyz      = ', opt_info%domain_xyz
      write(MsgOut,'(a)')     ''

    end if

    return

  end subroutine read_ctrl_boundary2

end module pr_control_mod
