!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   input_mod
!> @brief   input parameters and data for analysis
!! @authors Norio Takase (NT)
! 
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module input_mod

  use input_str_mod
  use fileio_control_mod
  use fileio_mt_mod
  use fileio_rtp_mod
  use fileio_atp_mod
  use fileio_grocrd_mod
  use fileio_grotop_mod
  use fileio_ambcrd_mod
  use fileio_prmtop_mod
  use fileio_top_mod
  use fileio_gpr_mod
  use fileio_pdb_mod
  use fileio_psf_mod
  use string_mod
  use messages_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_inp_info
    character(MaxFilename) :: ambcrdfile  = ''
    character(MaxFilename) :: ambreffile  = ''
    character(MaxFilename) :: atpfile     = ''
    character(MaxFilename) :: coorfile    = ''
    character(MaxFilename) :: cvfile      = ''
    character(MaxFilename) :: targetfile  = ''
    character(MaxFilename) :: dcdfile     = ''
    character(MaxFilename) :: dcdvelfile  = ''
    character(MaxFilename) :: excfile     = ''
    character(MaxFilename) :: gprfile     = ''
    character(MaxFilename) :: grocrdfile  = ''
    character(MaxFilename) :: groreffile  = ''
    character(MaxFilename) :: grotopfile  = ''
    character(MaxFilename) :: mtfile      = ''
    character(MaxFilename) :: indexfile   = ''
    character(MaxFilename) :: logfile     = ''
    character(MaxFilename) :: enefile     = ''
    character(MaxFilename) :: msdfile     = ''
    character(MaxFilename) :: pathfile    = ''
    character(MaxFilename) :: pathcvfile  = ''
    character(MaxFilename) :: pcafile     = ''
    character(MaxFilename) :: pdbfile     = ''
    character(MaxFilename) :: pdb_tgtfile = ''
    character(MaxFilename) :: pdb_avefile = ''
    character(MaxFilename) :: pdb_aftfile = ''
    character(MaxFilename) :: pdb_sphfile = ''
    character(MaxFilename) :: pdb_wbxfile = ''
    character(MaxFilename) :: prmtopfile  = ''
    character(MaxFilename) :: psffile     = ''
    character(MaxFilename) :: radfile     = ''
    character(MaxFilename) :: refenefile  = ''
    character(MaxFilename) :: reffile     = ''
    character(MaxFilename) :: fitfile     = ''
    character(MaxFilename) :: remfile     = ''
    character(MaxFilename) :: rstfile     = ''
    character(MaxFilename) :: rtpfile     = ''
    character(MaxFilename) :: topfile     = ''
    character(MaxFilename) :: valfile     = ''
    character(MaxFilename) :: vecfile     = ''
    character(MaxFilename) :: velfile     = ''
    character(MaxFilename) :: weightfile  = ''
    character(MaxFilename) :: xscfile     = ''
    character(MaxFilename) :: distfile    = ''
  end type s_inp_info

  public  :: show_ctrl_input
  public  :: read_ctrl_input
  public  :: setup_input
  public  :: input_files

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_input
  !> @brief        show control parameters in INPUT section
  !! @authors      NT
  !! @param[in]    keys : comma separated keys for show
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_input(keys)

    ! formal arguments
    character(*),            intent(in)    :: keys

    ! local variables
    integer                  :: i, ntok
    character(1000)          :: ckeys

    character(10), allocatable :: toks(:)


    if (len_trim(keys) == 0) &
      return

    ckeys = keys
    do i = 1, len(keys)
      if (ckeys(i:i) == ',') &
        ckeys(i:i) = ''
    end do

    ntok = split_num(ckeys)
    allocate(toks(ntok))
    call split(ntok, ntok, ckeys, toks)


    write(MsgOut,'(A)') '[INPUT]'

    do i = 1, ntok

      if (toks(i) == 'ambcrd') &
        write(MsgOut,'(A)') 'ambcrdfile     = input.crd       # AMBER coordinate file'

      if (toks(i) == 'ambref') &
        write(MsgOut,'(A)') 'ambreffile     = input.crd       # AMBER coordinate file'

      if (toks(i) == 'atp') &
        write(MsgOut,'(A)') 'atpfile        = input.atp       # GROMACS atp file'

      if (toks(i) == 'coor') &
        write(MsgOut,'(A)') 'coorfile       = input.coor      # coordinate restart file'

      if (toks(i) == 'cv') &
        write(MsgOut,'(A)') 'cvfile         = input.cv        # Collective variable file'

      if (toks(i) == 'target') &
        write(MsgOut,'(A)') 'targetfile     = input.target    # Target energy file'

      if (toks(i) == 'CV') &
        write(MsgOut,'(A)') 'cvfile         = input{}.cv      # Collective variable file'

      if (toks(i) == 'dcd') &
        write(MsgOut,'(A)') 'dcdfile        = input.dcd       # DCD file'

      if (toks(i) == 'DCD') &
        write(MsgOut,'(A)') 'dcdfile        = input{}.dcd     # DCD file'

      if (toks(i) == 'dcdvel') &
        write(MsgOut,'(A)') 'dcdvelfile     = input.dvl       # DCD velocity file'

      if (toks(i) == 'DCDVEL') &
        write(MsgOut,'(A)') 'dcdvelfile     = input{}.dvl     # DCD velocity file'

      if (toks(i) == 'exc') &
        write(MsgOut,'(A)') 'excfile        = input.exc       # macro-molecule exclusion list file'

!      if (toks(i) == 'gpr') &
!        write(MsgOut,'(A)') 'gprfile        = input.gpr       # GPR file'

      if (toks(i) == 'grocrd') &
        write(MsgOut,'(A)') 'grocrdfile     = input.crd       # GROMACS coordinate file'

      if (toks(i) == 'groref') &
        write(MsgOut,'(A)') 'groreffile     = input.crd       # GROMACS coordinate file'

      if (toks(i) == 'grotop') &
        write(MsgOut,'(A)') 'grotopfile     = input.top       # GROMACS topology file'

      if (toks(i) == 'mt') &
        write(MsgOut,'(A)') 'mtfile         = input.mt        # motion tree infom'

      if (toks(i) == 'index') &
        write(MsgOut,'(A)') 'indexfile      = input.idx       # Index file'

      if (toks(i) == 'LOG') &
        write(MsgOut,'(A)') 'logfile        = input{}.log     # REMD energy log file'

      if (toks(i) == 'ene') &
        write(MsgOut,'(A)') 'enefile        = input{}.ene     # gREST ene file'

      if (toks(i) == 'msd') &
        write(MsgOut,'(A)') 'msdfile        = input.msd       # Mean-square-displacement file'

      if (toks(i) == 'path') &
        write(MsgOut,'(A)') 'pathfile       = input.path      # PATH file'

      if (toks(i) == 'PATH') &
        write(MsgOut,'(A)') 'pathfile       = input{}.path    # PATH file'

      if (toks(i) == 'pathcv') &
        write(MsgOut,'(A)') 'pathcvfile     = input.pathcv    # PATH CV file'

      if (toks(i) == 'PATHCV') &
        write(MsgOut,'(A)') 'pathcvfile     = input{}.pathcv  # PATH CV file'

      if (toks(i) == 'pca') &
        write(MsgOut,'(A)') 'pcafile        = input.pca       # PCA file'

      if (toks(i) == 'pdb') &
        write(MsgOut,'(A)') 'pdbfile        = input.pdb       # PDB file'

      if (toks(i) == 'pdb_tgt') &
        write(MsgOut,'(A)') 'pdb_tgtfile    = input.pdb       # PDB file (Target coordinates for MultiBasin)'

      if (toks(i) == 'pdb_ave') &
        write(MsgOut,'(A)') 'pdb_avefile    = input_ave.pdb   # PDB file (Average coordinates)'

      if (toks(i) == 'pdb_aft') &
        write(MsgOut,'(A)') 'pdb_aftfile    = input_aft.pdb   # PDB file (Fitted Average coordinates)'

      if (toks(i) == 'pdb_sph') &
        write(MsgOut,'(A)') 'pdb_sphfile    = input.pdb       # PDB file (macro-molecule sphere coordinates)'

      if (toks(i) == 'pdb_wbx') &
        write(MsgOut,'(A)') 'pdb_wbxfile    = input.pdb       # PDB file (water box)'

      if (toks(i) == 'prmtop') &
        write(MsgOut,'(A)') 'prmtopfile     = input.top       # AMBER parameter topology file'

      if (toks(i) == 'psf') &
        write(MsgOut,'(A)') 'psffile        = input.psf       # protein structure file'

      if (toks(i) == 'rad') &
        write(MsgOut,'(A)') 'radfile        = input.rad       # radius file'

      if (toks(i) == 'ref') &
        write(MsgOut,'(A)') 'reffile        = input.pdb       # PDB file'

      if (toks(i) == 'fit') &
        write(MsgOut,'(A)') 'fitfile        = input.pdb       # PDB file'

      if (toks(i) == 'rem') &
        write(MsgOut,'(A)') 'remfile        = input.rem       # REMD parameter ID file'

      if (toks(i) == 'REM') &
        write(MsgOut,'(A)') 'remfile        = input{}.rem     # REMD parameter ID file'

      if (toks(i) == 'rst') &
        write(MsgOut,'(A)') 'rstfile        = input.rst       # restart file'

      if (toks(i) == 'RST') &
        write(MsgOut,'(A)') '# rstfile      = input{}.rst     # restart files'

      if (toks(i) == 'rtp') &
        write(MsgOut,'(A)') 'rtpfile        = input.rtp       # GROMACS rtp file'

      if (toks(i) == 'top') &
        write(MsgOut,'(A)') 'topfile        = input.top       # CHARMM topology file'

      if (toks(i) == 'val') &
        write(MsgOut,'(A)') 'valfile        = input.val       # VAL file'

      if (toks(i) == 'vec') &
        write(MsgOut,'(A)') 'vecfile        = input.vec       # VEC file'

      if (toks(i) == 'vel') &
        write(MsgOut,'(A)') 'velfile        = input.vel       # velocity restart file'

      if (toks(i) == 'weight') &
        write(MsgOut,'(A)') 'weightfile     = input.weight    # weight file'

      if (toks(i) == 'WEIGHT') &
        write(MsgOut,'(A)') 'weightfile     = input{}.weight  # weight file'

      if (toks(i) == 'xsc') &
        write(MsgOut,'(A)') 'xscfile        = input.xsc       # extended system file'

      if (toks(i) == 'dist') &
        write(MsgOut,'(A)') 'distfile       = input.dist       # pathcv distance file'

      if (toks(i) == 'DIST') &
        write(MsgOut,'(A)') 'distfile       = input{}.dist       # pathcv distance file'

    end do

    write(MsgOut,'(A)') ' '

    deallocate(toks)

    return

  end subroutine show_ctrl_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_input
  !> @brief        read control parameters in INPUT section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[out]   inp_info : INPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine read_ctrl_input(handle, inp_info)
  
    ! parameters
    character(*),            parameter :: Section = 'Input'
    
    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_inp_info),        intent(inout) :: inp_info


    ! read parameters from control file
    ! 
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'Ambcrdfile', &
         inp_info%ambcrdfile)
    call read_ctrlfile_string(handle, Section, 'Ambreffile', &
         inp_info%ambreffile)
    call read_ctrlfile_string(handle, Section, 'Atpfile', &
         inp_info%atpfile)
    call read_ctrlfile_string(handle, Section, 'Coorfile', &
         inp_info%coorfile)
    call read_ctrlfile_string(handle, Section, 'Cvfile', &
         inp_info%cvfile)
    call read_ctrlfile_string(handle, Section, 'Targetfile', &
         inp_info%targetfile)
    call read_ctrlfile_string(handle, Section, 'Dcdfile', &
         inp_info%dcdfile)
    call read_ctrlfile_string(handle, Section, 'Dcdvelfile', &
         inp_info%dcdvelfile)
    call read_ctrlfile_string(handle, Section, 'excfile', &
         inp_info%excfile)
    call read_ctrlfile_string(handle, Section, 'gprfile', &
         inp_info%gprfile)
    call read_ctrlfile_string(handle, Section, 'Groreffile', &
         inp_info%groreffile)
    call read_ctrlfile_string(handle, Section, 'Grocrdfile', &
         inp_info%grocrdfile)
    call read_ctrlfile_string(handle, Section, 'Groreffile', &
         inp_info%groreffile)
    call read_ctrlfile_string(handle, Section, 'Grotopfile', &
         inp_info%grotopfile)
    call read_ctrlfile_string(handle, Section, 'Mtfile', &
         inp_info%mtfile)
    call read_ctrlfile_string(handle, Section, 'Indexfile', &
         inp_info%indexfile)
    call read_ctrlfile_string(handle, Section, 'Logfile', &
         inp_info%logfile)
    call read_ctrlfile_string(handle, Section, 'Enefile', &
         inp_info%enefile)
    call read_ctrlfile_string(handle, Section, 'msdfile', &
         inp_info%msdfile)
    call read_ctrlfile_string(handle, Section, 'pathfile', &
         inp_info%pathfile)
    call read_ctrlfile_string(handle, Section, 'pathcvfile', &
         inp_info%pathcvfile)
    call read_ctrlfile_string(handle, Section, 'pcafile', &
         inp_info%pcafile)
    call read_ctrlfile_string(handle, Section, 'Pdbfile', &
         inp_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'Pdb_tgtfile', &
         inp_info%pdb_tgtfile)
    call read_ctrlfile_string(handle, Section, 'Pdb_avefile', &
         inp_info%pdb_avefile)
    call read_ctrlfile_string(handle, Section, 'Pdb_aftfile', &
         inp_info%pdb_aftfile)
    call read_ctrlfile_string(handle, Section, 'Pdb_sphfile', &
         inp_info%pdb_sphfile)
    call read_ctrlfile_string(handle, Section, 'Pdb_wbxfile', &
         inp_info%pdb_wbxfile)
    call read_ctrlfile_string(handle, Section, 'Prmtopfile', &
         inp_info%prmtopfile)
    call read_ctrlfile_string(handle, Section, 'Psffile', &
         inp_info%psffile)
    call read_ctrlfile_string(handle, Section, 'Radfile', &
         inp_info%radfile)
    call read_ctrlfile_string(handle, Section, 'Refenefile', &
         inp_info%refenefile)      
    call read_ctrlfile_string(handle, Section, 'Reffile', &
         inp_info%reffile)
    call read_ctrlfile_string(handle, Section, 'Fitfile', &
         inp_info%fitfile)
    call read_ctrlfile_string(handle, Section, 'Remfile', &
         inp_info%remfile)
    call read_ctrlfile_string(handle, Section, 'Rstfile', &
         inp_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'Rtpfile', &
         inp_info%rtpfile)
    call read_ctrlfile_string(handle, Section, 'Topfile', &
         inp_info%topfile)
    call read_ctrlfile_string(handle, Section, 'Valfile', &
         inp_info%valfile)
    call read_ctrlfile_string(handle, Section, 'Vecfile', &
         inp_info%vecfile)
    call read_ctrlfile_string(handle, Section, 'Velfile', &
         inp_info%velfile)
    call read_ctrlfile_string(handle, Section, 'Weightfile', &
         inp_info%weightfile)
    call read_ctrlfile_string(handle, Section, 'Xscfile', &
         inp_info%xscfile)
    call read_ctrlfile_string(handle, Section, 'Distfile', &
         inp_info%distfile)

    call end_ctrlfile_section(handle)


    ! write parameters to MsgOut
    !
    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Input> Input Files'
      if (len_trim(inp_info%ambcrdfile) > 0) &
        write(MsgOut,'(A20,A)') '  ambcrdfile      = ', trim(inp_info%ambcrdfile)
      if (len_trim(inp_info%ambreffile) > 0) &
        write(MsgOut,'(A20,A)') '  ambreffile      = ', trim(inp_info%ambreffile)
      if (len_trim(inp_info%atpfile) > 0) &
        write(MsgOut,'(A20,A)') '  atpfile         = ', trim(inp_info%atpfile)
      if (len_trim(inp_info%coorfile) > 0) &
        write(MsgOut,'(A20,A)') '  coorfile        = ', trim(inp_info%coorfile)
      if (len_trim(inp_info%cvfile) > 0) &
        write(MsgOut,'(A20,A)') '  cvfile          = ', trim(inp_info%cvfile)
      if (len_trim(inp_info%targetfile) > 0) &
        write(MsgOut,'(A20,A)') '  targetfile      = ', trim(inp_info%targetfile)
      if (len_trim(inp_info%dcdfile) > 0) &
        write(MsgOut,'(A20,A)') '  dcdfile         = ', trim(inp_info%dcdfile)
      if (len_trim(inp_info%dcdvelfile) > 0) &
        write(MsgOut,'(A20,A)') '  dcdvelfile      = ', trim(inp_info%dcdvelfile)
      if (len_trim(inp_info%excfile) > 0) &
        write(MsgOut,'(A20,A)') '  excfile         = ', trim(inp_info%excfile)
      if (len_trim(inp_info%gprfile) > 0) &
        write(MsgOut,'(A20,A)') '  gprfile         = ', trim(inp_info%gprfile)
      if (len_trim(inp_info%groreffile) > 0) &
        write(MsgOut,'(A20,A)') '  groreffile      = ', trim(inp_info%groreffile)
      if (len_trim(inp_info%grocrdfile) > 0) &
        write(MsgOut,'(A20,A)') '  grocrdfile      = ', trim(inp_info%grocrdfile)
      if (len_trim(inp_info%grotopfile) > 0) &
        write(MsgOut,'(A20,A)') '  grotopfile      = ', trim(inp_info%grotopfile)
      if (len_trim(inp_info%indexfile) > 0) &
        write(MsgOut,'(A20,A)') '  indexfile       = ', trim(inp_info%indexfile)
      if (len_trim(inp_info%logfile) > 0) &
        write(MsgOut,'(A20,A)') '  logfile         = ', trim(inp_info%logfile)
      if (len_trim(inp_info%enefile) > 0) &
        write(MsgOut,'(A20,A)') '  enefile         = ', trim(inp_info%enefile)
      if (len_trim(inp_info%msdfile) > 0) &
        write(MsgOut,'(A20,A)') '  msdfile         = ', trim(inp_info%msdfile)
      if (len_trim(inp_info%mtfile) > 0) &
        write(MsgOut,'(A20,A)') '  mtfile          = ', trim(inp_info%mtfile)
      if (len_trim(inp_info%pathfile) > 0) &
        write(MsgOut,'(A20,A)') '  pathfile        = ', trim(inp_info%pathfile)
      if (len_trim(inp_info%pathcvfile) > 0) &
        write(MsgOut,'(A20,A)') '  pathcvfile      = ', trim(inp_info%pathcvfile)
      if (len_trim(inp_info%pcafile) > 0) &
        write(MsgOut,'(A20,A)') '  pcafile         = ', trim(inp_info%pcafile)
      if (len_trim(inp_info%pdbfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdbfile         = ', trim(inp_info%pdbfile)
      if (len_trim(inp_info%pdb_avefile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_avefile     = ', trim(inp_info%pdb_avefile)
      if (len_trim(inp_info%pdb_aftfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_aftfile     = ', trim(inp_info%pdb_aftfile)
      if (len_trim(inp_info%pdb_sphfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_sphfile     = ', trim(inp_info%pdb_sphfile)
      if (len_trim(inp_info%pdb_wbxfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_wbxfile     = ', trim(inp_info%pdb_wbxfile)
      if (len_trim(inp_info%prmtopfile) > 0) &
        write(MsgOut,'(A20,A)') '  prmtopfile      = ', trim(inp_info%prmtopfile)
      if (len_trim(inp_info%psffile) > 0) &
        write(MsgOut,'(A20,A)') '  psffile         = ', trim(inp_info%psffile)
      if (len_trim(inp_info%radfile) > 0) &
        write(MsgOut,'(A20,A)') '  radfile         = ', trim(inp_info%radfile)
      if (len_trim(inp_info%refenefile) > 0) &
        write(MsgOut,'(A20,A)') '  refenefile      = ', trim(inp_info%refenefile)        
      if (len_trim(inp_info%reffile) > 0) &
        write(MsgOut,'(A20,A)') '  reffile         = ', trim(inp_info%reffile)
      if (len_trim(inp_info%fitfile) > 0) &
        write(MsgOut,'(A20,A)') '  fitfile         = ', trim(inp_info%fitfile)
      if (len_trim(inp_info%remfile) > 0) &
        write(MsgOut,'(A20,A)') '  remfile         = ', trim(inp_info%remfile)
      if (len_trim(inp_info%rstfile) > 0) &
        write(MsgOut,'(A20,A)') '  rstfile         = ', trim(inp_info%rstfile)
      if (len_trim(inp_info%rtpfile) > 0) &
        write(MsgOut,'(A20,A)') '  rtpfile         = ', trim(inp_info%rtpfile)
      if (len_trim(inp_info%topfile) > 0) &
        write(MsgOut,'(A20,A)') '  topfile         = ', trim(inp_info%topfile)
      if (len_trim(inp_info%valfile) > 0) &
        write(MsgOut,'(A20,A)') '  valfile         = ', trim(inp_info%valfile)
      if (len_trim(inp_info%vecfile) > 0) &
        write(MsgOut,'(A20,A)') '  vecfile         = ', trim(inp_info%vecfile)
      if (len_trim(inp_info%velfile) > 0) &
        write(MsgOut,'(A20,A)') '  velfile         = ', trim(inp_info%velfile)
      if (len_trim(inp_info%weightfile) > 0) &
        write(MsgOut,'(A20,A)') '  weightfile      = ', trim(inp_info%weightfile)
      if (len_trim(inp_info%xscfile) > 0) &
        write(MsgOut,'(A20,A)') '  xscfile         = ', trim(inp_info%xscfile)
      if (len_trim(inp_info%distfile) > 0) &
        write(MsgOut,'(A20,A)') '  distfile        = ', trim(inp_info%distfile)
      write(MsgOut,'(A)') ' '
    end if

    return

  end subroutine read_ctrl_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_input
  !> @brief        setup input files
  !! @authors      NT
  !! @param[in]    inp_info : INPUT section control parameters information
  !! @param[out]   input    : input information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_input(inp_info, input)

    ! formal argments
    type(s_inp_info),        intent(in)    :: inp_info
    type(s_input),           intent(inout) :: input


    input%ambcrdfile  = inp_info%ambcrdfile
    input%ambreffile  = inp_info%ambreffile
    input%atpfile     = inp_info%atpfile
    input%coorfile    = inp_info%coorfile
    input%cvfile      = inp_info%cvfile
    input%targetfile  = inp_info%targetfile
    input%dcdfile     = inp_info%dcdfile
    input%dcdvelfile  = inp_info%dcdvelfile
    input%excfile     = inp_info%excfile
    input%gprfile     = inp_info%gprfile
    input%grocrdfile  = inp_info%grocrdfile
    input%groreffile  = inp_info%groreffile
    input%grotopfile  = inp_info%grotopfile
    input%mtfile      = inp_info%mtfile
    input%indexfile   = inp_info%indexfile
    input%logfile     = inp_info%logfile
    input%enefile     = inp_info%enefile
    input%msdfile     = inp_info%msdfile
    input%pathfile    = inp_info%pathfile
    input%pathcvfile  = inp_info%pathcvfile
    input%pcafile     = inp_info%pcafile
    input%pdbfile     = inp_info%pdbfile
    input%pdb_tgtfile = inp_info%pdb_tgtfile
    input%pdb_avefile = inp_info%pdb_avefile
    input%pdb_aftfile = inp_info%pdb_aftfile
    input%pdb_sphfile = inp_info%pdb_sphfile
    input%pdb_wbxfile = inp_info%pdb_wbxfile
    input%prmtopfile  = inp_info%prmtopfile
    input%psffile     = inp_info%psffile
    input%radfile     = inp_info%radfile
    input%refenefile  = inp_info%refenefile    
    input%reffile     = inp_info%reffile
    input%fitfile     = inp_info%fitfile
    input%remfile     = inp_info%remfile
    input%rstfile     = inp_info%rstfile
    input%rtpfile     = inp_info%rtpfile
    input%topfile     = inp_info%topfile
    input%valfile     = inp_info%valfile
    input%vecfile     = inp_info%vecfile
    input%velfile     = inp_info%velfile
    input%weightfile  = inp_info%weightfile
    input%xscfile     = inp_info%xscfile
    input%distfile    = inp_info%distfile

    return

  end subroutine setup_input

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    input_files
  !> @brief        read input data from files
  !! @authors      NT
  !! @param[in]    inp_info : INPUT section control parameters information
  !! @param[out]   psf      : PSF information [optional]
  !! @param[out]   ref      : coordinate data for reference [optional]
  !! @param[out]   pdb      : PDB information [optional]
  !! @param[out]   gpr      : GPR information [optional]
  !! @param[out]   top      : CHARMM topology information [optional]
  !! @param[out]   prmtop   : AMBER parameter / topology information [optional]
  !! @param[out]   ambcrd   : AMBER coordinate information [optional]
  !! @param[out]   ambref   : AMBER coordinate information [optional]
  !! @param[out]   grotop   : GROMACS parameter / topology information [optiona]
  !! @param[out]   grocrd   : GROMACS coordinate information [optional]
  !! @param[out]   groref   : GROMACS coordinate information [optional]
  !! @param[out]   atp      : GROMACS atp information [optional]
  !! @param[out]   rtp      : GROMACS rtp information [optional]
  !! @param[out]   mt       : Motion Tree information [optional]
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine input_files(inp_info, psf, ref, fit, pdb, gpr, top, prmtop, &
                         ambref, ambcrd, grotop, groref, grocrd, atp, rtp, mt)

    ! formal arguments
    type(s_inp_info),         intent(in)    :: inp_info
    type(s_psf),    optional, intent(inout) :: psf
    type(s_pdb),    optional, intent(inout) :: ref
    type(s_pdb),    optional, intent(inout) :: fit
    type(s_pdb),    optional, intent(inout) :: pdb
    type(s_gpr),    optional, intent(inout) :: gpr
    type(s_top),    optional, intent(inout) :: top
    type(s_prmtop), optional, intent(inout) :: prmtop
    type(s_ambcrd), optional, intent(inout) :: ambcrd
    type(s_ambcrd), optional, intent(inout) :: ambref
    type(s_grotop), optional, intent(inout) :: grotop
    type(s_grocrd), optional, intent(inout) :: groref
    type(s_grocrd), optional, intent(inout) :: grocrd
    type(s_atp),    optional, intent(inout) :: atp
    type(s_rtp),    optional, intent(inout) :: rtp
    type(s_mt),     optional, intent(inout) :: mt


    if (present(gpr))    call init_gpr(gpr)
    if (present(psf))    call init_psf(psf)
    if (present(prmtop)) call init_prmtop(prmtop)
    if (present(top))    call init_top(top)
    if (present(grotop)) call init_grotop(grotop)
    if (present(pdb))    call init_pdb(pdb)
    if (present(fit))    call init_pdb(fit)
    if (present(ambcrd)) call init_ambcrd(ambcrd)
    if (present(ambref)) call init_ambcrd(ambref)
    if (present(grocrd)) call init_grocrd(grocrd)
    if (present(groref)) call init_grocrd(groref)
    if (present(ref))    call init_pdb(ref)

    if (present(gpr))    call init_gpr(gpr)
    if (present(psf))    call init_psf(psf)
    if (present(prmtop)) call init_prmtop(prmtop)
    if (present(grotop)) call init_grotop(grotop)
    if (present(pdb))    call init_pdb(pdb)
    if (present(fit))    call init_pdb(fit)
    if (present(ambcrd)) call init_ambcrd(ambcrd)
    if (present(grocrd)) call init_grocrd(grocrd)
    if (present(ref))    call init_pdb(ref)

    if (len_trim(inp_info%psffile) > 0 .and. present(psf)) then
      call input_psf(inp_info%psffile, psf)
    end if

    if (len_trim(inp_info%reffile) > 0 .and. present(ref)) then
      call input_pdb(inp_info%reffile, ref)
    end if

    if (len_trim(inp_info%fitfile) > 0 .and. present(fit)) then
      call input_pdb(inp_info%fitfile, fit)
    end if

    if (len_trim(inp_info%pdbfile) > 0 .and. present(pdb)) then
      call input_pdb(inp_info%pdbfile, pdb)
    end if

    if (len_trim(inp_info%gprfile) > 0 .and. present(gpr)) then
      call input_gpr(inp_info%gprfile, gpr)
    end if

    if (len_trim(inp_info%prmtopfile) > 0 .and. present(prmtop)) then
      call input_prmtop(inp_info%prmtopfile, prmtop)
    end if

    if (len_trim(inp_info%ambreffile) > 0 .and. present(ambref)) then
      call input_ambcrd(inp_info%ambreffile, ambcrd)
    end if

    if (len_trim(inp_info%ambcrdfile) > 0 .and. present(ambcrd)) then
      call input_ambcrd(inp_info%ambcrdfile, ambcrd)
    end if

    if (len_trim(inp_info%ambreffile) > 0 .and. present(ambref)) then
      call input_ambcrd(inp_info%ambreffile, ambref)
    end if

    if (len_trim(inp_info%grotopfile) > 0 .and. present(grotop)) then
      call input_grotop(inp_info%grotopfile, grotop)
    end if

    if (len_trim(inp_info%groreffile) > 0 .and. present(groref)) then
      call input_grocrd(inp_info%groreffile, grocrd)
    end if

    if (len_trim(inp_info%grocrdfile) > 0 .and. present(grocrd)) then
      call input_grocrd(inp_info%grocrdfile, grocrd)
    end if

    if (len_trim(inp_info%groreffile) > 0 .and. present(groref)) then
      call input_grocrd(inp_info%groreffile, groref)
    end if

    if (len_trim(inp_info%atpfile) > 0 .and. present(atp)) then
      call input_atp(inp_info%atpfile, atp)
    end if

    if (len_trim(inp_info%mtfile) > 0 .and. present(mt)) then
      call input_mt(inp_info%mtfile, mt)
    end if

    if (len_trim(inp_info%rtpfile) > 0 .and. present(rtp)) then
      call input_rtp(inp_info%rtpfile, rtp)
    end if

    if (len_trim(inp_info%topfile) > 0 .and. present(top)) then
      call input_top(inp_info%topfile, top)
    end if

    return

  end subroutine input_files

end module input_mod
