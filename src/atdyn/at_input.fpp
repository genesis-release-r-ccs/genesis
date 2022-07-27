!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   at_input_mod
!> @brief   read parameters and data for md simulation
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Chigusa Kobayashi (CK), 
!!          Jaewoon Jung (JJ), Norio Takase (NT), Kiyoshi Yagi (KY)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module at_input_mod

  use at_remd_str_mod
  use at_rpath_str_mod
  use fileio_table_mod
  use fileio_mode_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_psf_mod
  use fileio_rst_mod
  use fileio_gpr_mod
  use fileio_str_mod
  use fileio_par_mod
  use fileio_top_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_mode_mod
  use fileio_eef1_mod
  use fileio_morph_mod
  use fileio_spot_mod
  use fileio_control_mod
  use fileio_rstmep_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_inp_info
    character(MaxFilenameLong) :: topfile    = ''
    character(MaxFilenameLong) :: parfile    = ''
    character(MaxFilenameLong) :: strfile    = ''
    character(MaxFilename)     :: gprfile    = ''
    character(MaxFilename)     :: psffile    = ''
    character(MaxFilename)     :: prmtopfile = ''
    character(MaxFilename)     :: grotopfile = ''
    character(MaxFilename)     :: pdbfile    = ''
    character(MaxFilename)     :: crdfile    = ''
    character(MaxFilename)     :: ambcrdfile = ''
    character(MaxFilename)     :: grocrdfile = ''
    character(MaxFilename)     :: rstfile    = ''
    character(MaxFilename)     :: reffile    = ''
    character(MaxFilename)     :: fitfile    = ''
    character(MaxFilename)     :: ambreffile = ''
    character(MaxFilename)     :: groreffile = ''
    character(MaxFilename)     :: modefile   = ''
    character(MaxFilename)     :: rstmepfile = ''
    character(MaxFilename)     :: eef1file   = ''
    character(MaxFilename)     :: tablefile  = ''
    character(MaxFilename)     :: morphfile  = ''
    character(MaxFilename)     :: spotfile   = ''
    character(MaxFilename)     :: minfofile   = ''
  end type s_inp_info

  ! subroutines
  public  :: show_ctrl_input
  public  :: read_ctrl_input
  public  :: input_md
  public  :: input_min
  public  :: input_remd
  public  :: input_rpath
  public  :: input_bd
  public  :: input_morph
  public  :: input_rpath_resetup
  private :: include_id_to_filename

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_input
  !> @brief        show INPUT section usage
  !! @authors      NT
  !! @param[in]    show_all : show all usage or not
  !! @param[in]    run_mode : run mode string : "md", "min", "remd", "rpath"
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_input(show_all, run_mode)

    ! formal arguments
    logical,                 intent(in)    :: show_all
    character(*),            intent(in)    :: run_mode


    if (show_all) then

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  CHARMM Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample.crd      # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample.rst      # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# eef1file = sample.inp     # solvation parameter file for EEF1/IMM1'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '# spotfile = sample.xyz     # center of spherical potential'
        write(MsgOut,'(A)') '##  Go-model'
        write(MsgOut,'(A)') '# psffile = sample.psf      # CHARMM protein structure file'
        write(MsgOut,'(A)') '# gprfile = sample.gpr      # GENESIS Go-model parameter file'
        write(MsgOut,'(A)') '# pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample.crd      # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample.rst      # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference PDB for positional restraint'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        write(MsgOut,'(A)') '##  Lookup table file'
        write(MsgOut,'(A)') '# tablefile = sample.tbl    # lookup table file'
        write(MsgOut,'(A)') ' '


      case ('min', 'vib')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  CHARMM Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample.crd      # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample.rst      # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# eef1file = sample.inp     # solvation parameter file for EEF1/IMM1'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '# spotfile = sample.xyz     # center of spherical potential'
        write(MsgOut,'(A)') '##  Go-model'
        write(MsgOut,'(A)') '# psffile = sample.psf      # CHARMM protein structure file'
        write(MsgOut,'(A)') '# gprfile = sample.gpr      # GENESIS Go-model parameter file'
        write(MsgOut,'(A)') '# pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample.crd      # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample.rst      # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference PDB for positional restraint'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        write(MsgOut,'(A)') '##  Minfo file for restart'
        write(MsgOut,'(A)') '# minfofile = sample.minfo  # file obtained from a less number of vibatom'
        write(MsgOut,'(A)') ' '


      case ('remd')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  CHARMM Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample{}.pdb    # PDB file'
        write(MsgOut,'(A)') '# crdfile  = sample{}.crd    # coordinates file'
        write(MsgOut,'(A)') '# rstfile  = sample{}.rst    # restart file'
        write(MsgOut,'(A)') '# reffile  = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# eef1file = sample.inp     # solvation parameter file for EEF1/IMM1'
        write(MsgOut,'(A)') '##  Lookup table file'
        write(MsgOut,'(A)') '# tablefile = sample.tbl    # lookup table file'
        write(MsgOut,'(A)') '# spotfile = sample.xyz     # center of spherical potential'
        write(MsgOut,'(A)') ' '


      case ('rpath')
        
        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '## Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample{}.pdb    # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample{}.crd    # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample{}.rst    # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# eef1file = sample.inp     # solvation parameter file for EEF1/IMM1'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '# spotfile = sample.xyz     # center of spherical potential'
        write(MsgOut,'(A)') '##  Go-model'
        write(MsgOut,'(A)') '# psffile = sample.psf      # CHARMM protein structure file'
        write(MsgOut,'(A)') '# gprfile = sample.gpr      # GENESIS Go-model parameter file'
        write(MsgOut,'(A)') '# pdbfile = sample{}.pdb    # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample{}.crd    # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample{}.rst    # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference PDB for positional restraint'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        !write(MsgOut,'(A)') '##  Rpath (MEP, FEP)'
        !write(MsgOut,'(A)') '# rstmepfile = sample{}.rstmep  # restart file for MEP and FEP'
        write(MsgOut,'(A)') ' '

        write(MsgOut,'(A)') '##  Lookup table file'
        write(MsgOut,'(A)') '# tablefile = sample.tbl    # lookup table file'
        write(MsgOut,'(A)') ' '

      case ('bd')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') 'grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') 'grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        write(MsgOut,'(A)') '##  Lookup table file'
        write(MsgOut,'(A)') '# tablefile = sample.tbl    # lookup table file'
        write(MsgOut,'(A)') 'morphfile = sample.mor     # morph file'
        write(MsgOut,'(A)') ' '

      end select

    else

      select case (run_mode)

      case ('md')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') ' '


      case ('min', 'vib')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') ' '


      case ('remd')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') 'topfile  =  sample.top     # topology file'
        write(MsgOut,'(A)') 'parfile  =  sample.par     # parameter file'
        write(MsgOut,'(A)') 'psffile  =  sample.psf     # protein structure file'
        write(MsgOut,'(A)') 'pdbfile  =  sample{}.pdb   # PDB file'
        write(MsgOut,'(A)') ' '


      case ('rpath')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') 'topfile  =  sample.top     # topology file'
        write(MsgOut,'(A)') 'parfile  =  sample.par     # parameter file'
        write(MsgOut,'(A)') 'psffile  =  sample.psf     # protein structure file'
        write(MsgOut,'(A)') 'pdbfile  =  sample{}.pdb   # PDB file'
        write(MsgOut,'(A)') ' '

      case ('bd')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') 'grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') 'grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# morphfile = sample.mor      # Morph file'
        write(MsgOut,'(A)') ' '

      end select

    end if

    return

  end subroutine show_ctrl_input
  
  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_input
  !> @brief        read INPUT section in the control file
  !! @authors      YS, TM, JJ, CK
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   inp_info : INPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_input(handle, inp_info)
  
    ! parameters
    character(*),            parameter     :: Section = 'Input'

    ! formal arguments
    integer,                 intent(in)    :: handle
    type(s_inp_info),        intent(inout) :: inp_info


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'Topfile',   inp_info%topfile)
    call read_ctrlfile_string(handle, Section, 'Parfile',   inp_info%parfile)
    call read_ctrlfile_string(handle, Section, 'Strfile',   inp_info%strfile)
    call read_ctrlfile_string(handle, Section, 'Gprfile',   inp_info%gprfile)
    call read_ctrlfile_string(handle, Section, 'Psffile',   inp_info%psffile)
    call read_ctrlfile_string(handle, Section, 'Prmtopfile',inp_info%prmtopfile)
    call read_ctrlfile_string(handle, Section, 'Grotopfile',inp_info%grotopfile)
    call read_ctrlfile_string(handle, Section, 'Pdbfile',   inp_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'Crdfile',   inp_info%crdfile)
    call read_ctrlfile_string(handle, Section, 'Ambcrdfile',inp_info%ambcrdfile)
    call read_ctrlfile_string(handle, Section, 'Grocrdfile',inp_info%grocrdfile)
    call read_ctrlfile_string(handle, Section, 'Rstfile',   inp_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'Reffile',   inp_info%reffile)
    call read_ctrlfile_string(handle, Section, 'Fitfile',   inp_info%fitfile)
    call read_ctrlfile_string(handle, Section, 'Ambreffile',inp_info%ambreffile)
    call read_ctrlfile_string(handle, Section, 'Groreffile',inp_info%groreffile)
    call read_ctrlfile_string(handle, Section, 'Modefile',  inp_info%modefile)
    call read_ctrlfile_string(handle, Section, 'rstmepfile',inp_info%rstmepfile)
    call read_ctrlfile_string(handle, Section, 'Eef1file',  inp_info%eef1file)
    call read_ctrlfile_string(handle, Section, 'Tablefile', inp_info%tablefile)
    call read_ctrlfile_string(handle, Section, 'Morphfile', inp_info%morphfile)
    call read_ctrlfile_string(handle, Section, 'Spotfile',  inp_info%spotfile)
    call read_ctrlfile_string(handle, Section, 'minfofile', inp_info%minfofile)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Input> Input Files'

      if (inp_info%topfile .ne. '') then
        write(MsgOut,*) ' topfile = ', trim(inp_info%topfile)
      end if

      if (inp_info%parfile .ne. '') then
        write(MsgOut,*) ' parfile = ', trim(inp_info%parfile)
      end if

      if (inp_info%strfile .ne. '') then
        write(MsgOut,*) ' strfile = ', trim(inp_info%strfile)
      end if

      if (inp_info%gprfile .ne. '') then
        write(MsgOut,*) ' gprfile = ', trim(inp_info%gprfile)
      end if

      if (inp_info%psffile .ne. '') then
        write(MsgOut,*) ' psffile = ', trim(inp_info%psffile)
      end if

      if (inp_info%prmtopfile .ne. '') then
        write(MsgOut,*) ' prmtopfile = ', trim(inp_info%prmtopfile)
      end if

      if (inp_info%grotopfile .ne. '') then
        write(MsgOut,*) ' grotopfile = ', trim(inp_info%grotopfile)
      end if

      if (inp_info%pdbfile .ne. '') then
        write(MsgOut,*) ' pdbfile = ', trim(inp_info%pdbfile)
      end if

      if (inp_info%crdfile .ne. '') then
        write(MsgOut,*) ' crdfile = ', trim(inp_info%crdfile)
      end if

      if (inp_info%ambcrdfile .ne. '') then
        write(MsgOut,*) ' ambcrdfile = ', trim(inp_info%ambcrdfile)
      end if

      if (inp_info%grocrdfile .ne. '') then
        write(MsgOut,*) ' grocrdfile = ', trim(inp_info%grocrdfile)
      end if

      if (inp_info%rstfile .ne. '') then
        write(MsgOut,*) ' rstfile = ', trim(inp_info%rstfile)
      end if

      if (inp_info%reffile .ne. '') then
        write(MsgOut,*) ' reffile = ', trim(inp_info%reffile)
      end if

      if (inp_info%fitfile .ne. '') then
        write(MsgOut,*) ' fitfile = ', trim(inp_info%fitfile)
      end if

      if (inp_info%ambreffile .ne. '') then
        write(MsgOut,*) ' ambreffile = ', trim(inp_info%ambreffile)
      end if

      if (inp_info%groreffile .ne. '') then
        write(MsgOut,*) ' groreffile = ', trim(inp_info%groreffile)
      end if

      if (inp_info%modefile .ne. '') then
        write(MsgOut,*) ' modefile = ', trim(inp_info%modefile)
      end if

      if (inp_info%eef1file .ne. '') then
        write(MsgOut,*) ' eef1file = ', trim(inp_info%eef1file)
      end if

      if (inp_info%rstmepfile .ne. '') then
        write(MsgOut,*) ' rstmepfile = ', trim(inp_info%rstmepfile)
      end if

      if (inp_info%tablefile .ne. '') then
        write(MsgOut,*) ' tablefile = ', trim(inp_info%tablefile)
      end if

      if (inp_info%spotfile .ne. '') then
        write(MsgOut,*) ' spotfile = ', trim(inp_info%spotfile)
      end if

      if (inp_info%morphfile .ne. '') then
        write(MsgOut,*) ' morphfile = ', trim(inp_info%morphfile)
      end if

      if (inp_info%minfofile .ne. '') then
        write(MsgOut,*) ' minfofile = ', trim(inp_info%minfofile)
      end if

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_md
  !> @brief        read input data from files
  !! @authors      YS, TM, CK, NT
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
  !! @param[out]   gpr      : information of GENESIS Go model parameters
  !! @param[out]   psf      : information of Protein structure data
  !! @param[out]   prmtop   : information of AMBER parameter topology data
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode     : information of principal component vector
  !! @param[out]   table    : information of lookup table data
  !! @param[out]   spot     : information of spherical potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_md(inp_info, top, par, gpr, psf, prmtop, grotop,       &
                      pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                      mode, eef1, table, morph_in, spot)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_eef1),            intent(inout) :: eef1
    type(s_table),           intent(inout) :: table
    type(s_morph_in),        intent(inout) :: morph_in
    type(s_spot),            intent(inout) :: spot


    call init_top(top)
    call init_par(par)
    call init_gpr(gpr)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_rst(rst)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)
    call init_morph_in(morph_in)

    if (inp_info%topfile .ne. '') then
      call input_top(inp_info%topfile, top)
    end if

    if (inp_info%parfile .ne. '') then
      call input_par(inp_info%parfile, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      call input_str(inp_info%strfile, top, par)
    end if

    if (inp_info%gprfile .ne. '') then
      call input_gpr(inp_info%gprfile, gpr)
    end if

    if (inp_info%psffile .ne. '') then
      call input_psf(inp_info%psffile, psf)
    end if

    if (inp_info%prmtopfile .ne. '') then
      call input_prmtop(inp_info%prmtopfile, prmtop)
    end if

    if (inp_info%grotopfile .ne. '') then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (inp_info%pdbfile .ne. '') then
      call input_pdb(inp_info%pdbfile, pdb)
    end if

    if (inp_info%crdfile .ne. '') then
      call input_crd(inp_info%crdfile, crd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      call input_ambcrd(inp_info%ambcrdfile, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (inp_info%rstfile .ne. '') then
      call input_rst(inp_info%rstfile, rst)
    end if

    if (inp_info%reffile .ne. '') then
      call input_pdb(inp_info%reffile, ref)
    end if

    if (inp_info%ambreffile .ne. '') then
      call input_ambcrd(inp_info%ambreffile, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (inp_info%modefile .ne. '') then
      call input_mode(inp_info%modefile, mode)
    end if

    if (inp_info%eef1file .ne. '') then
      call input_eef1(inp_info%eef1file, eef1)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    if (inp_info%morphfile .ne. '') then
      call input_morph_in(inp_info%morphfile, morph_in)
    end if

    if (inp_info%spotfile .ne. '') then
      call input_spot(inp_info%spotfile, spot)
    end if

    return

  end subroutine input_md

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_min
  !> @brief        read input data from files
  !! @authors      TM, CK
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
  !! @param[out]   gpr      : information of GENESIS Go model parameters
  !! @param[out]   psf      : information of Protein structure data
  !! @param[out]   prmtop   : information of AMBER parameter topology data
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode     : information of principal component vector
  !! @param[out]   table    : information of lookup table data
  !! @param[out]   spot     : information of spherical potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_min(inp_info, top, par, gpr, psf, prmtop, grotop,       &
                       pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                       mode, eef1, table, morph_in, spot)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_eef1),            intent(inout) :: eef1
    type(s_table),           intent(inout) :: table
    type(s_morph_in),        intent(inout) :: morph_in
    type(s_spot),            intent(inout) :: spot


    call init_top(top)
    call init_par(par)
    call init_gpr(gpr)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_rst(rst)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)
    call init_morph_in(morph_in)

    if (inp_info%topfile .ne. '') then
      call input_top(inp_info%topfile, top)
    end if

    if (inp_info%parfile .ne. '') then
      call input_par(inp_info%parfile, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      call input_str(inp_info%strfile, top, par)
    end if

    if (inp_info%gprfile .ne. '') then
      call input_gpr(inp_info%gprfile, gpr)
    end if

    if (inp_info%psffile .ne. '') then
      call input_psf(inp_info%psffile, psf)
    end if

    if (inp_info%prmtopfile .ne. '') then
      call input_prmtop(inp_info%prmtopfile, prmtop)
    end if

    if (inp_info%grotopfile .ne. '') then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (inp_info%pdbfile .ne. '') then
      call input_pdb(inp_info%pdbfile, pdb)
    end if

    if (inp_info%crdfile .ne. '') then
      call input_crd(inp_info%crdfile, crd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      call input_ambcrd(inp_info%ambcrdfile, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (inp_info%rstfile .ne. '') then
      call input_rst(inp_info%rstfile, rst)
    end if

    if (inp_info%reffile .ne. '') then
      call input_pdb(inp_info%reffile, ref)
    end if

    if (inp_info%ambreffile .ne. '') then
      call input_ambcrd(inp_info%ambreffile, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (inp_info%modefile .ne. '') then
      call input_mode(inp_info%modefile, mode)
    end if

    if (inp_info%eef1file .ne. '') then
      call input_eef1(inp_info%eef1file, eef1)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    if (inp_info%morphfile .ne. '') then
      call input_morph_in(inp_info%morphfile, morph_in)
    end if

    if (inp_info%spotfile .ne. '') then
      call input_spot(inp_info%spotfile, spot)
    end if

    return

  end subroutine input_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_remd
  !> @brief        read input data from files
  !! @authors      TM, CK
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
  !! @param[out]   gpr      : information of GENESIS Go model parameters
  !! @param[out]   psf      : information of Protein structure data
  !! @param[out]   prmtop   : information of AMBER parameter topology data
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode     : information of principal component vector
  !! @param[out]   table    : information of lookup table data
  !! @param[out]   spot     : information of spherical potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_remd(inp_info, top, par, gpr, psf, prmtop, grotop, pdb, &
                        crd, ambcrd, grocrd, rst, ref, ambref, groref,     &
                        mode, eef1, table, morph_in, spot)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_eef1),            intent(inout) :: eef1
    type(s_table),           intent(inout) :: table
    type(s_morph_in),        intent(inout) :: morph_in
    type(s_spot),            intent(inout) :: spot

    ! local variables
    character(MaxMultiFilename)            :: filename
    character(MaxFilenameLong)             :: filename_long


    call init_top(top)
    call init_par(par)
    call init_gpr(gpr)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_rst(rst)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)
    call init_morph_in(morph_in)

    if (inp_info%topfile .ne. '') then
      filename_long = inp_info%topfile
      call input_top(filename_long, top)
    end if

    if (inp_info%parfile .ne. '') then
      filename_long = inp_info%parfile
      call input_par(filename_long, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      filename_long = inp_info%strfile
      call input_str(filename_long, top, par)
    end if

    if (inp_info%gprfile .ne. '') then
      filename = inp_info%gprfile
      call input_gpr(filename, gpr)
    end if

    if (inp_info%psffile .ne. '') then
      filename = inp_info%psffile
      call input_psf(filename, psf)
    end if

    if (inp_info%prmtopfile .ne. '') then
      call input_prmtop(inp_info%prmtopfile, prmtop)
    end if

    if (inp_info%grotopfile .ne. '') then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (inp_info%reffile .ne. '') then
      filename = inp_info%reffile
      call include_id_to_filename(filename)
      call input_pdb(filename, ref)
    end if

    if (inp_info%ambreffile .ne. '') then
      filename = inp_info%ambreffile
      call input_ambcrd(filename, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      filename = inp_info%groreffile
      call input_grocrd(filename, groref)
    end if

    if (inp_info%crdfile .ne. '') then
      filename = inp_info%crdfile
      call include_id_to_filename(filename)
      call input_crd(filename, crd)
    end if

    if (inp_info%pdbfile .ne. '') then
      filename = inp_info%pdbfile
      call include_id_to_filename(filename)
      call input_pdb(filename, pdb)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      filename = inp_info%ambcrdfile
      call include_id_to_filename(filename)
      call input_ambcrd(filename, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      filename = inp_info%grocrdfile
      call include_id_to_filename(filename)
      call input_grocrd(filename, grocrd)
    end if

    rst%rstfile_type = RstfileTypeUndef
    if (inp_info%rstfile .ne. '') then
      filename = inp_info%rstfile
      call include_id_to_filename(filename)
      call input_rst(filename, rst)
    end if

    if (inp_info%modefile .ne. '') then
      filename = inp_info%modefile
      call input_mode(filename, mode)
    end if

    if (inp_info%eef1file .ne. '') then
      call input_eef1(inp_info%eef1file, eef1)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    if (inp_info%morphfile .ne. '') then
      call input_morph_in(inp_info%morphfile, morph_in)
    end if

    if (inp_info%spotfile .ne. '') then
      call input_spot(inp_info%spotfile, spot)
    end if

    return

  end subroutine input_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rpath
  !> @brief        read input data from files
  !! @authors      TM, CK, YK, YM, KY
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
  !! @param[out]   gpr      : information of GENESIS Go model parameters
  !! @param[out]   psf      : information of Protein structure data
  !! @param[out]   prmtop   : information of AMBER parameter topology data
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   fit      : information of fit coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode     : information of principal component vector
  !! @param[out]   rstmep   : information of rpath
  !! @param[out]   table    : information of lookup table data
  !! @param[out]   spot     : information of spherical potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rpath(inp_info, top, par, gpr, psf, prmtop, grotop, pdb, &
                         crd, ambcrd, grocrd, rst, ref, fit, ambref,        &
                         groref, mode, rstmep, eef1, table, morph_in, spot)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_pdb),             intent(inout) :: fit
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_rstmep),          intent(inout) :: rstmep
    type(s_eef1),            intent(inout) :: eef1
    type(s_table),           intent(inout) :: table
    type(s_morph_in),        intent(inout) :: morph_in
    type(s_spot),            intent(inout) :: spot

    ! local variables
    character(MaxMultiFilename)            :: filename
    character(MaxFilenameLong)             :: filename_long


    call init_top(top)
    call init_par(par)
    call init_gpr(gpr)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_rst(rst)
    call init_pdb(ref)
    call init_pdb(fit)
    call init_ambcrd(ambref)
    call init_grocrd(groref)
    call init_morph_in(morph_in)

    if (inp_info%topfile .ne. '') then
      filename_long = inp_info%topfile
      call input_top(filename_long, top)
    end if

    if (inp_info%parfile .ne. '') then
      filename_long = inp_info%parfile
      call input_par(filename_long, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      filename_long = inp_info%strfile
      call input_str(filename_long, top, par)
    end if

    if (inp_info%gprfile .ne. '') then
      filename = inp_info%gprfile
      call input_gpr(filename, gpr)
    end if

    if (inp_info%psffile .ne. '') then
      filename = inp_info%psffile
      call input_psf(filename, psf)
    end if

    if (inp_info%prmtopfile .ne. '') then
      filename = inp_info%prmtopfile
      call input_prmtop(filename, prmtop)
    end if

    if (inp_info%grotopfile .ne. '') then
      filename = inp_info%grotopfile
      call input_grotop(filename, grotop)
    end if

    if (inp_info%pdbfile .ne. '') then
      filename = inp_info%pdbfile
      call include_id_to_filename(filename)
      call input_pdb(filename, pdb)
    end if

    if (inp_info%crdfile .ne. '') then
      filename = inp_info%crdfile
      call include_id_to_filename(filename)
      call input_crd(filename, crd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      filename = inp_info%ambcrdfile
      call include_id_to_filename(filename)
      call input_ambcrd(filename, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      filename = inp_info%grocrdfile
      call include_id_to_filename(filename)
      call input_grocrd(filename, grocrd)
    end if

    rst%rstfile_type = RstfileTypeUndef
    if (inp_info%rstfile .ne. '') then
      filename = inp_info%rstfile
      call include_id_to_filename(filename)
      call input_rst(filename, rst)
    end if

    if (inp_info%reffile .ne. '') then
      filename = inp_info%reffile
      call include_id_to_filename(filename)
      call input_pdb(filename, ref)
    end if

    if (inp_info%rstmepfile .ne. '') then
      filename = inp_info%rstmepfile
      call include_id_to_filename(filename)
      call input_rstmep(filename, rstmep)
    end if

    if (inp_info%fitfile .ne. '') then
      filename = inp_info%fitfile
      call input_pdb(filename, fit)
    end if

    if (inp_info%ambreffile .ne. '') then
      filename = inp_info%ambreffile
      call include_id_to_filename(filename)
      call input_ambcrd(filename, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      filename = inp_info%groreffile
      call include_id_to_filename(filename)
      call input_grocrd(filename, groref)
    end if

    if (inp_info%modefile .ne. '') then
      filename = inp_info%modefile
      call include_id_to_filename(filename)
      call input_mode(filename, mode)
    end if


    if (inp_info%eef1file .ne. '') then
      call input_eef1(inp_info%eef1file, eef1)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    if (inp_info%morphfile .ne. '') then
      call input_morph_in(inp_info%morphfile, morph_in)
    end if

    if (inp_info%spotfile .ne. '') then
      call input_spot(inp_info%spotfile, spot)
    end if

    return

  end subroutine input_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rpath_resetup
  !> @brief        read input data from files for the case "nreplica > nproc"
  !! @authors      SI
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   fit      : information of fit coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode     : information of principal component vector
  !! @param[out]   rstmep   : information of rpath
  !! @param[out]   spot     : information of spherical potential
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rpath_resetup(inp_info, top, par, gpr, psf, &
                                 prmtop, grotop, pdb, crd, ambcrd, grocrd,    &
                                 rst, ref, fit, ambref, groref, mode, rstmep, &
                                 eef1, spot, rpath, rst_filename)

    ! formal arguments
    type(s_inp_info),        intent(inout) :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_pdb),             intent(inout) :: fit
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_rstmep),          intent(inout) :: rstmep
    type(s_eef1),            intent(inout) :: eef1
    type(s_spot),            intent(inout) :: spot
    type(s_rpath),           intent(inout) :: rpath

    character(MaxFilename),  intent(in)    :: rst_filename
    ! local variables
    character(MaxMultiFilename)            :: filename
    character(MaxFilenameLong)             :: filename_long


    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (rpath%first_replica .and. rpath%first_iter) then

      call init_top(top)
      call init_par(par)
      call init_gpr(gpr)
      call init_psf(psf)
      call init_pdb(fit)
      call init_prmtop(prmtop)
      call init_grotop(grotop)

      if (inp_info%topfile .ne. '') then
        filename_long = inp_info%topfile
        call input_top(filename_long, top)
      end if

      if (inp_info%parfile .ne. '') then
        filename_long = inp_info%parfile
        call input_par(filename_long, par, top)
      end if

      if (inp_info%strfile .ne. '') then
        filename_long = inp_info%strfile
        call input_str(filename_long, top, par)
      end if

      if (inp_info%gprfile .ne. '') then
        filename = inp_info%gprfile
        call input_gpr(filename, gpr)
      end if

      if (inp_info%psffile .ne. '') then
        filename = inp_info%psffile
        call input_psf(filename, psf)
      end if

      if (inp_info%prmtopfile .ne. '') then
        filename = inp_info%prmtopfile
        call input_prmtop(filename, prmtop)
      end if

      if (inp_info%grotopfile .ne. '') then
        filename = inp_info%grotopfile
        call input_grotop(filename, grotop)
      end if

      if (inp_info%fitfile .ne. '') then
        filename = inp_info%fitfile
        call input_pdb(filename, fit)
      end if

      if (inp_info%eef1file .ne. '') then
        call input_eef1(inp_info%eef1file, eef1)
      end if

      if (inp_info%spotfile .ne. '') then
        call input_spot(inp_info%spotfile, spot)
      end if
      !rpath%first_replica = .false.
    end if

    if (inp_info%pdbfile .ne. '') then
      filename = inp_info%pdbfile
      call include_id_to_filename(filename)
      call input_pdb(filename, pdb)
    end if

    if (inp_info%crdfile .ne. '') then
      filename = inp_info%crdfile
      call include_id_to_filename(filename)
      call input_crd(filename, crd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      filename = inp_info%ambcrdfile
      call include_id_to_filename(filename)
      call input_ambcrd(filename, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      filename = inp_info%grocrdfile
      call include_id_to_filename(filename)
      call input_grocrd(filename, grocrd)
    end if

    rst%rstfile_type = RstfileTypeUndef
    if (.not. rpath%first_iter) inp_info%rstfile = rst_filename
    if (inp_info%rstfile .ne. '') then
      filename = inp_info%rstfile
      call include_id_to_filename(filename)
      call input_rst(filename, rst)
    end if

    if (inp_info%reffile .ne. '') then
      filename = inp_info%reffile
      call include_id_to_filename(filename)
      call input_pdb(filename, ref)
    end if

    if (inp_info%rstmepfile .ne. '') then
      filename = inp_info%rstmepfile
      call include_id_to_filename(filename)
      call input_rstmep(filename, rstmep)
    end if

    if (inp_info%ambreffile .ne. '') then
      filename = inp_info%ambreffile
      call include_id_to_filename(filename)
      call input_ambcrd(filename, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      filename = inp_info%groreffile
      call include_id_to_filename(filename)
      call input_grocrd(filename, groref)
    end if

    if (inp_info%modefile .ne. '') then
      filename = inp_info%modefile
      call include_id_to_filename(filename)
      call input_mode(filename, mode)
    end if

    return

  end subroutine input_rpath_resetup

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_bd
  !> @brief        read input data from files
  !! @authors      TA
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   table    : information of lookup table data
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_bd(inp_info, grotop, grocrd, groref, rst, table)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_grotop),          intent(inout) :: grotop
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_grocrd),          intent(inout) :: groref
    type(s_rst),             intent(inout) :: rst
    type(s_table),           intent(inout) :: table


    call init_grotop(grotop)
    call init_grocrd(grocrd)
    call init_grocrd(groref)
    call init_rst(rst)

    if (inp_info%grotopfile .ne. '') then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (inp_info%grocrdfile .ne. '') then
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (inp_info%groreffile .ne. '') then
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (inp_info%rstfile .ne. '') then
      call input_rst(inp_info%rstfile, rst)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    return

  end subroutine input_bd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_morph
  !> @brief        read input data from files
  !! @authors      TM, CK
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
  !! @param[out]   gpr      : information of GENESIS Go model parameters
  !! @param[out]   psf      : information of Protein structure data
  !! @param[out]   prmtop   : information of AMBER parameter topology data
  !! @param[out]   grotop   : information of GROMACS parameter topology data
  !! @param[out]   pdb      : information of coordinate data
  !! @param[out]   crd      : information of coordinate data
  !! @param[out]   ambcrd   : information of AMBER coordinate data
  !! @param[out]   grocrd   : information of GROMACS coordinate data
  !! @param[out]   rst      : information of restart data
  !! @param[out]   ref      : information of reference coordinate data
  !! @param[out]   ambref   : information of refernece AMBER coordinate data
  !! @param[out]   groref   : information of refernece GROMACS coordinate data
  !! @param[out]   mode 
  !! @param[out]   morph    : information of morphing file
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_morph(inp_info, top, par, gpr, psf, prmtop, grotop,        &
                         pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref,  &
                         mode, eef1, table, morph_in)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
    type(s_gpr),             intent(inout) :: gpr
    type(s_psf),             intent(inout) :: psf
    type(s_prmtop),          intent(inout) :: prmtop
    type(s_grotop),          intent(inout) :: grotop
    type(s_pdb),             intent(inout) :: pdb
    type(s_crd),             intent(inout) :: crd
    type(s_ambcrd),          intent(inout) :: ambcrd
    type(s_grocrd),          intent(inout) :: grocrd
    type(s_rst),             intent(inout) :: rst
    type(s_pdb),             intent(inout) :: ref
    type(s_ambcrd),          intent(inout) :: ambref
    type(s_grocrd),          intent(inout) :: groref
    type(s_mode),            intent(inout) :: mode
    type(s_eef1),            intent(inout) :: eef1
    type(s_table),           intent(inout) :: table
    type(s_morph_in),        intent(inout) :: morph_in


    call init_top(top)
    call init_par(par)
    call init_gpr(gpr)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_rst(rst)
    call init_ambcrd(ambref)
    call init_grocrd(groref)
    call init_morph_in(morph_in)

    if (inp_info%topfile .ne. '') then
      call input_top(inp_info%topfile, top)
    end if

    if (inp_info%parfile .ne. '') then
      call input_par(inp_info%parfile, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      call input_str(inp_info%strfile, top, par)
    end if

    if (inp_info%gprfile .ne. '') then
      call input_gpr(inp_info%gprfile, gpr)
    end if

    if (inp_info%psffile .ne. '') then
      call input_psf(inp_info%psffile, psf)
    end if

    if (inp_info%prmtopfile .ne. '') then
      call input_prmtop(inp_info%prmtopfile, prmtop)
    end if

    if (inp_info%grotopfile .ne. '') then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (inp_info%pdbfile .ne. '') then
      call input_pdb(inp_info%pdbfile, pdb)
    end if

    if (inp_info%crdfile .ne. '') then
      call input_crd(inp_info%crdfile, crd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      call input_ambcrd(inp_info%ambcrdfile, ambcrd)
    end if

    if (inp_info%grocrdfile .ne. '') then
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (inp_info%rstfile .ne. '') then
      call input_rst(inp_info%rstfile, rst)
    end if

    if (inp_info%reffile .ne. '') then
      call input_pdb(inp_info%reffile, ref)
    end if

    if (inp_info%ambreffile .ne. '') then
      call input_ambcrd(inp_info%ambreffile, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (inp_info%modefile .ne. '') then
      call input_mode(inp_info%modefile, mode)
    end if

    if (inp_info%eef1file .ne. '') then
      call input_eef1(inp_info%eef1file, eef1)
    end if

    if (inp_info%tablefile .ne. '') then
      call input_table(inp_info%tablefile, table)
    end if

    if (inp_info%morphfile .ne. '') then
      call input_morph_in(inp_info%morphfile, morph_in)
    end if

    return

  end subroutine input_morph

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    include_id_to_filename
  !> @brief        include id to filename
  !! @authors      TM
  !! @param[in]    id       : index
  !! @param[in]    ndigit   : number of digit
  !! @param[inout] filename : replicate filename
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine include_id_to_filename(filename)

    ! formal arguments
    character(MaxFilename),  intent(inout) :: filename

    ! local variables
    integer                  :: comp, ndigit, id
    integer                  :: i, j, ci1, ci2, cnumber
    character(MaxFilename)   :: filename_ori
    character(10)            :: frmt, cd, cid


    ! define replica id
    !
    id = my_country_no + 1
    if (nrep_per_proc > 1) id = my_replica_no
    do i = 1, 100
      comp = 10**i
      if (id < comp) then
        ndigit = i
        exit
      end if
    end do

    ! check filename
    !
    filename_ori = filename

    ci1 = scan(filename, '{')
    ci2 = scan(filename, '}')

    if (ci1 == 0 .or. ci2 == 0) &
      return

    if (ci1 > 0 .and. ci2 > ci1) then

      write(cd,'(i10)') ndigit
      frmt = '(i' // trim(adjustl(cd)) // '.' // trim(adjustl(cd)) // ')'
      write(cid,frmt) id

      cnumber = len_trim(filename_ori)
      if (cnumber + ndigit > MaxFilename) &
         call error_msg('Error: too long filename'//filename_ori)

      j = 0
      do i = 1, ci1 - 1
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do
      do i = 1, ndigit
        j = j + 1
        filename(j:j) = cid(i:i)
      end do
      do i = ci2+1, MaxFilename
        j = j + 1
        filename(j:j) = filename_ori(i:i)
      end do

    end if

    return

  end subroutine include_id_to_filename

end module at_input_mod
