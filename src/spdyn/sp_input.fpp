!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_input_mod
!> @brief   read parameters and data for md simulation
!! @authors Yuji Sugita (YS), Takaharu Mori (TM), Jaewoon Jung (JJ), 
!!          Chigusa Kobayashi (CK), Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_input_mod

  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_localres_mod
  use fileio_psf_mod
  use fileio_rst_mod
  use fileio_str_mod
  use fileio_par_mod
  use fileio_top_mod
  use fileio_crd_mod
  use fileio_pdb_mod
  use fileio_mode_mod
  use fileio_control_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_inp_info
    character(MaxFilename) :: topfile    = ''
    character(MaxFilename) :: parfile    = ''
    character(MaxFilenameLong) :: strfile    = ''
    character(MaxFilename) :: psffile    = ''
    character(MaxFilename) :: prmtopfile = ''
    character(MaxFilename) :: grotopfile = ''
    character(MaxFilename) :: pdbfile    = ''
    character(MaxFilename) :: crdfile    = ''
    character(MaxFilename) :: ambcrdfile = ''
    character(MaxFilename) :: grocrdfile = ''
    character(MaxFilename) :: rstfile    = ''
    character(MaxFilename) :: selfile    = ''
    character(MaxFilename) :: reffile    = ''
    character(MaxFilename) :: fitfile    = ''
    character(MaxFilename) :: ambreffile = ''
    character(MaxFilename) :: groreffile = ''
    character(MaxFilename) :: local_resfile = ''
    character(MaxFilename) :: modefile   = ''
  end type s_inp_info

  ! subroutines
  public  :: show_ctrl_input
  public  :: read_ctrl_input
  public  :: input_md
  public  :: input_min
  public  :: input_remd
  public  :: input_rpath
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
        write(MsgOut,'(A)') '# selfile = sample.sel      # selection list file (used only for parallel I/O'
        write(MsgOut,'(A)') '# rstfile = sample().rst    # domain restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '# localresfile = sample.txt # local restraint file'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        write(MsgOut,'(A)') ' '


      case ('min')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  CHARMM Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'pdbfile = sample.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample.crd      # CHARMM coordinates file'
        write(MsgOut,'(A)') '# selfile = sample.sel      # selction list file'
        write(MsgOut,'(A)') '# rstfile = sample.rst      # GENESIS restart file'
        write(MsgOut,'(A)') '# rstfile = sample().rst    # domain restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '# localresfile = sample.txt # local restraint file'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
        write(MsgOut,'(A)') ' '


      case ('remd')

        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '##  CHARMM Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'selfile = sample.sel      # selection list file'
        write(MsgOut,'(A)') 'pdbfile = sample{}.pdb      # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample{}.crd    # coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample{}.rst    # restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# localresfile = sample.txt # local restraint file'
        write(MsgOut,'(A)') ' '


      case ('rpath')
        
        write(MsgOut,'(A)') '[INPUT]'
        write(MsgOut,'(A)') '## Force Field'
        write(MsgOut,'(A)') 'topfile = sample.top      # topology file'
        write(MsgOut,'(A)') 'parfile = sample.par      # parameter file'
        write(MsgOut,'(A)') 'strfile = sample.str      # stream file'
        write(MsgOut,'(A)') 'psffile = sample.psf      # protein structure file'
        write(MsgOut,'(A)') 'selfile = sample.sel      # selection list file'
        write(MsgOut,'(A)') 'pdbfile = sample{}.pdb    # PDB file'
        write(MsgOut,'(A)') '# crdfile = sample{}.crd    # CHARMM coordinates file'
        write(MsgOut,'(A)') '# rstfile = sample{}.rst    # GENESIS restart file'
        write(MsgOut,'(A)') '# reffile = sample.pdb      # reference for restraint'
        write(MsgOut,'(A)') '# modefile = sample.mode    # principal component vector'
        write(MsgOut,'(A)') '##  AMBER Force Field'
        write(MsgOut,'(A)') '# prmtopfile = sample.top   # AMBER parameter topology file'
        write(MsgOut,'(A)') '# ambcrdfile = sample.crd   # AMBER coordinate file'
        write(MsgOut,'(A)') '# ambreffile = sample.crd   # reference AMBER coordinate file'
        write(MsgOut,'(A)') '##  GROMACS Force Field'
        write(MsgOut,'(A)') '# grotopfile = sample.top   # GROMACS parameter topology file'
        write(MsgOut,'(A)') '# grocrdfile = sample.crd   # GROMACS coordinate file'
        write(MsgOut,'(A)') '# groreffile = sample.crd   # reference GROMACS coordinate file'
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


      case ('min')

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
    call read_ctrlfile_string(handle, Section, 'Psffile',   inp_info%psffile)
    call read_ctrlfile_string(handle, Section, 'Prmtopfile',inp_info%prmtopfile)
    call read_ctrlfile_string(handle, Section, 'Grotopfile',inp_info%grotopfile)
    call read_ctrlfile_string(handle, Section, 'Pdbfile',   inp_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'Crdfile',   inp_info%crdfile)
    call read_ctrlfile_string(handle, Section, 'Ambcrdfile',inp_info%ambcrdfile)
    call read_ctrlfile_string(handle, Section, 'Grocrdfile',inp_info%grocrdfile)
    call read_ctrlfile_string(handle, Section, 'selfile',   inp_info%selfile)
    call read_ctrlfile_string(handle, Section, 'Rstfile',   inp_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'Reffile',   inp_info%reffile)
    call read_ctrlfile_string(handle, Section, 'Fitfile',   inp_info%fitfile)
    call read_ctrlfile_string(handle, Section, 'Ambreffile',inp_info%ambreffile)
    call read_ctrlfile_string(handle, Section, 'Groreffile',inp_info%groreffile)
    call read_ctrlfile_string(handle, Section, 'Modefile',  inp_info%modefile)
    call read_ctrlfile_string(handle, Section, 'Localresfile', &
                              inp_info%local_resfile)

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
        write(MsgOut,*) ' selfile = ', trim(inp_info%selfile)
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

      if (inp_info%local_resfile .ne. '') then
        write(MsgOut,*) ' localresfile = ', trim(inp_info%local_resfile)
      end if

      if (inp_info%modefile .ne. '') then
        write(MsgOut,*) ' modefile = ', trim(inp_info%modefile)
      end if

      write(MsgOut,'(A)') ' '

    end if

    return

  end subroutine read_ctrl_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_md
  !> @brief        read input data from files
  !! @authors      YS, TM, CK
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
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
  !! @param[out]   localres : information of local restraint data
  !! @param[out]   mode     : information of principal component vector
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_md(inp_info, top, par, psf, prmtop, grotop, pdb, crd,  &
                      ambcrd, grocrd, rst, ref, ambref, groref, localres, &
                      mode)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
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
    type(s_localres),        intent(inout) :: localres
    type(s_mode),            intent(inout) :: mode


    call init_top(top)
    call init_par(par)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (inp_info%topfile .ne. '') then
      call input_top(inp_info%topfile, top)
    end if

    if (inp_info%parfile .ne. '') then
      call input_par(inp_info%parfile, par)
    end if

    if (inp_info%strfile .ne. '') then
      call input_str(inp_info%strfile, top, par)
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

    if (inp_info%local_resfile .ne. '') then
      call input_localres(inp_info%local_resfile, localres)
    end if

    if (inp_info%modefile .ne. '') then
      call input_mode(inp_info%modefile, mode)
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
  !! @param[out]   localres : information of local restraint data
  !! @param[out]   mode     : information of principal component vector
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_min(inp_info, top, par, psf, prmtop, grotop, pdb, crd,  &
                       ambcrd, grocrd, rst, ref, ambref, groref, localres, &
                       mode)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
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
    type(s_localres),        intent(inout) :: localres
    type(s_mode),            intent(inout) :: mode


    call init_top(top)
    call init_par(par)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (inp_info%topfile .ne. '') then
      call input_top(inp_info%topfile, top)
    end if

    if (inp_info%parfile .ne. '') then
      call input_par(inp_info%parfile, par)
    end if

    if (inp_info%strfile .ne. '') then
      call input_str(inp_info%strfile, top, par)
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

    if (inp_info%local_resfile .ne. '') then
      call input_localres(inp_info%local_resfile, localres)
    end if

    if (inp_info%modefile .ne. '') then
      call input_mode(inp_info%modefile, mode)
    end if

    return

  end subroutine input_min

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_remd
  !> @brief        read input data from files
  !! @authors      YS, TM, CK
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
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
  !! @param[out]   localres : information of local restraint data
  !! @param[out]   mode     : information of principal component vector
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_remd(inp_info, top, par, psf, prmtop, grotop, &
                        pdb, crd, ambcrd, grocrd, rst, ref, ambref, groref, &
                        localres, mode)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
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
    type(s_localres),        intent(inout) :: localres
    type(s_mode),            intent(inout) :: mode

    ! local variables
    character(MaxMultiFilename)            :: filename
    character(MaxFilenameLong)             :: filename_long


    call init_top(top)
    call init_par(par)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (inp_info%topfile .ne. '') then
      filename = inp_info%topfile
      call input_top(filename, top)
    end if

    if (inp_info%parfile .ne. '') then
      filename = inp_info%parfile
      call input_par(filename, par)
    end if

    if (inp_info%strfile .ne. '') then
      filename_long = inp_info%strfile
      call input_str(filename_long, top, par)
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

    if (inp_info%grocrdfile .ne. '') then
      filename = inp_info%grocrdfile
      call include_id_to_filename(filename)
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (inp_info%ambcrdfile .ne. '') then
      filename = inp_info%ambcrdfile
      call include_id_to_filename(filename)
      call input_ambcrd(inp_info%ambcrdfile, ambcrd)
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

    if (inp_info%ambreffile .ne. '') then
      filename = inp_info%ambreffile
      call include_id_to_filename(filename)
      call input_ambcrd(inp_info%ambreffile, ambref)
    end if

    if (inp_info%groreffile .ne. '') then
      filename = inp_info%groreffile
      call include_id_to_filename(filename)
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (inp_info%local_resfile .ne. '') then
      filename = inp_info%local_resfile
      call input_localres(filename, localres)
    end if

    if (inp_info%modefile .ne. '') then
      filename = inp_info%modefile
      call input_mode(filename, mode)
    end if

    return

  end subroutine input_remd

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_rpath
  !> @brief        read input data from files
  !! @authors      TM, CK, YK, YM
  !! @param[in]    inp_info : information of INPUT section control parameters
  !! @param[out]   top      : information of CHARMM topology data
  !! @param[out]   par      : information of CHARMM force field parameters
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
  !! @param[out]   localres : information of local restraint data
  !! @param[out]   mode     : information of principal component vector
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_rpath(inp_info, top, par, psf, prmtop, grotop, pdb,  &
                         crd, ambcrd, grocrd, rst, ref, fit, ambref, groref, &
                         localres, mode)

    ! formal arguments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_top),             intent(inout) :: top
    type(s_par),             intent(inout) :: par
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
    type(s_localres),        intent(inout) :: localres
    type(s_mode),            intent(inout) :: mode

    ! local variables
    character(MaxMultiFilename)            :: filename
    character(MaxFilenameLong)             :: filename_long


    call init_top(top)
    call init_par(par)
    call init_psf(psf)
    call init_prmtop(prmtop)
    call init_grotop(grotop)
    call init_pdb(pdb)
    call init_crd(crd)
    call init_ambcrd(ambcrd)
    call init_grocrd(grocrd)
    call init_pdb(ref)
    call init_pdb(fit)
    call init_ambcrd(ambref)
    call init_grocrd(groref)

    if (inp_info%topfile .ne. '') then
      filename = inp_info%topfile
      call input_top(filename, top)
    end if

    if (inp_info%parfile .ne. '') then
      filename = inp_info%parfile
      call input_par(filename, par, top)
    end if

    if (inp_info%strfile .ne. '') then
      filename_long = inp_info%strfile
      call input_str(filename_long, top, par)
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

    if (inp_info%local_resfile .ne. '') then
      filename = inp_info%local_resfile
      call include_id_to_filename(filename)
      call input_localres(filename, localres)
    end if

    if (inp_info%modefile .ne. '') then
      filename = inp_info%modefile
      call include_id_to_filename(filename)
      call input_mode(filename, mode)
    end if

    return

  end subroutine input_rpath

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    include_id_to_filename
  !> @brief        include id to filename
  !! @authors      TM
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

    if (ci1 == 0 .or. ci2 ==0) &
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

end module sp_input_mod
