!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   output_mod
!> @brief   output data for analysis
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module output_mod

  use constants_mod
  use output_str_mod
  use fileio_control_mod
  use fileio_mod
  use string_mod
  use messages_mod
  use getopt_mod
  use mpi_parallel_mod

  implicit none
  private

  ! structures
  type, public :: s_out_info
    character(MaxFilename) :: ambcrdfile   = ''
    character(MaxFilename) :: angfile      = ''
    character(MaxFilename) :: cntfile      = ''
    character(MaxFilename) :: comangfile   = ''
    character(MaxFilename) :: comdisfile   = ''
    character(MaxFilename) :: comtorfile   = ''
    character(MaxFilename) :: coorfile     = ''
    character(MaxFilename) :: crdfile      = ''
    character(MaxFilename) :: crsfile      = ''
    character(MaxFilename) :: disfile      = ''
    character(MaxFilename) :: enefile      = ''
    character(MaxFilename) :: excfile      = ''
    character(MaxFilename) :: fenefile     = ''
    character(MaxFilename) :: gprfile      = ''
    character(MaxFilename) :: grotopfile   = ''
    character(MaxFilename) :: grocrdfile   = ''
    character(MaxFilename) :: grocrd_tgtfile   = ''
    character(MaxFilename) :: hb_listfile  = ''
    character(MaxFilename) :: indexfile    = ''
    character(MaxFilename) :: logfile      = ''
    character(MaxFilename) :: mapfile      = ''
    character(MaxFilename) :: morphfile    = ''
    character(MaxFilename) :: msdfile      = ''
    character(MaxFilename) :: outfile      = ''
    character(MaxFilename) :: parfile      = ''
    character(MaxFilename) :: pathcvfile   = ''
    character(MaxFilename) :: pcafile      = ''
    character(MaxFilename) :: pdbfile      = ''
    character(MaxFilename) :: pdb_avefile  = ''
    character(MaxFilename) :: pdb_aftfile  = ''
    character(MaxFilename) :: pdb_tgtfile  = ''
    character(MaxFilename) :: pmffile      = ''
    character(MaxFilename) :: pmlfile      = ''
    character(MaxFilename) :: prjfile      = ''
    character(MaxFilename) :: probfile     = ''
    character(MaxFilename) :: qmmm_crdfile = ''
    character(MaxFilename) :: qmmm_psffile = ''
    character(MaxFilename) :: qmmm_pdbfile = ''
    character(MaxFilename) :: qntfile      = ''
    character(MaxFilename) :: rdffile      = ''
    character(MaxFilename) :: rmsfile      = ''
    character(MaxFilename) :: rgfile       = ''
    character(MaxFilename) :: rstfile      = ''
    character(MaxFilename) :: tblfile      = ''
    character(MaxFilename) :: txtfile      = ''
    character(MaxFilename) :: topfile      = ''
    character(MaxFilename) :: torfile      = ''
    character(MaxFilename) :: trjfile      = ''
    character(MaxFilename) :: trrfile      = ''
    character(MaxFilename) :: valfile      = ''
    character(MaxFilename) :: vcvfile      = ''
    character(MaxFilename) :: vecfile      = ''
    character(MaxFilename) :: velfile      = ''
    character(MaxFilename) :: voronoifile  = ''
    character(MaxFilename) :: vmdfile      = ''
    character(MaxFilename) :: weightfile   = ''
    character(MaxFilename) :: xscfile      = ''
!TM(181230)
!    integer                :: allow_backup = tristate_NOT_SET
!TM(181230)
  end type s_out_info

  ! subroutines
  public  :: show_ctrl_output
  public  :: read_ctrl_output
  public  :: parse_ctrl_output
  public  :: setup_output

contains

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    show_ctrl_output
  !> @brief        show control parameters in OUTPUT section
  !! @authors      NT
  !! @param[in]    keys : comma separated keys for show
  !
  !======1=========2=========3=========4=========5=========6=========7=========8
  
  subroutine show_ctrl_output(keys)

    ! formal arguments
    character(*),            intent(in)    :: keys

    ! local variables
    integer                  :: i, ntok
    character(Maxline)       :: ckeys
!TM(181230)
!    character(3)             :: yes_no
!TM(181230)
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


    write(MsgOut,'(A)') '[OUTPUT]'

    do i = 1, ntok

      if (toks(i) == 'ambcrd') &
        write(MsgOut,'(A)') 'ambcrdfile     = output.inpcrd    # AMBER CRD file'

      if (toks(i) == 'AMBCRD') &
        write(MsgOut,'(A)') '# ambcrdfile   = output.inpcrd.{} # AMBER CRD files'

      if (toks(i) == 'ang') &
        write(MsgOut,'(A)') 'angfile        = output.ang      # angle file'

      if (toks(i) == 'cnt') &
        write(MsgOut,'(A)') 'cntfile        = output.cnt      # CNT file'

      if (toks(i) == 'comang') &
        write(MsgOut,'(A)') 'comangfile     = output.comang   # COM angle file'

      if (toks(i) == 'comdis') &
        write(MsgOut,'(A)') 'comdisfile     = output.comdis   # COM distance file'

      if (toks(i) == 'comtor') &
        write(MsgOut,'(A)') 'comtorfile     = output.comtor   # COM torsion file'

      if (toks(i) == 'coor') &
        write(MsgOut,'(A)') 'coorfile       = output.coor     # coordinate restart file'

      if (toks(i) == 'crd') &
        write(MsgOut,'(A)') 'crdfile        = output.crd      # CHARMM CRD file'

      if (toks(i) == 'CRD') &
        write(MsgOut,'(A)') '# crdfile      = output.crd_{}   # CHARMM CRD files'

      if (toks(i) == 'crs') &
        write(MsgOut,'(A)') 'crsfile        = output.crs      # CRS file'

      if (toks(i) == 'dis') &
        write(MsgOut,'(A)') 'disfile        = output.dis      # distance file'

      if (toks(i) == 'ene') &
        write(MsgOut,'(A)') 'enefile        = output.ene      # REMD energy file'

      if (toks(i) == 'exc') &
        write(MsgOut,'(A)') 'excfile        = output.exc      # macro-molecule exclusion list file'

      if (toks(i) == 'fene') &
        write(MsgOut,'(A)') 'fenefile       = output.fene     # free energy file'

      if (toks(i) == 'map') &
       write(MsgOut,'(A)')  'mapfile        = output.xplor    # 3D density map file (xplor/ccp4/dx/sit)'

      if (toks(i) == 'msd') &
       write(MsgOut,'(A)')  'msdfile        = output.msd      # mean square displacement file'

      if (toks(i) == 'grotop') &
        write(MsgOut,'(A)') 'grotopfile     = output.grotop   # GROMACS topology file'

      if (toks(i) == 'grocrd') &
        write(MsgOut,'(A)') 'grocrdfile     = output.grocrd   # GROMACS coordinate file'

      if (toks(i) == 'index') &
        write(MsgOut,'(A)') 'indexfile      = output.idx      # Index file'

      if (toks(i) == 'LOG') &
        write(MsgOut,'(A)') 'logfile        = output{}.log    # REMD energy log file'

      if (toks(i) == 'out') &
        write(MsgOut,'(A)') 'outfile        = output.out      # output file'

      if (toks(i) == 'par') &
        write(MsgOut,'(A)') 'parfile        = output.par      # parameter file'

      if (toks(i) == 'grocrd_tgt') &
        write(MsgOut,'(A)') 'grocrd_tgtfile  = output.grocrd   # GROMACS coordinate file'

      if (toks(i) == 'pathcv') &
        write(MsgOut,'(A)') 'pathcvfile     = output.pathcv   # PATH CV file'

      if (toks(i) == 'PATHCV') &
        write(MsgOut,'(A)') 'pathcvfile     = output{}.pathcv # PATH CV file'

      if (toks(i) == 'pca') &
        write(MsgOut,'(A)') 'pcafile        = output.pca      # PCA file'

      if (toks(i) == 'pdb') &
        write(MsgOut,'(A)') 'pdbfile        = output.pdb      # PDB file'

      if (toks(i) == 'PDB') &
        write(MsgOut,'(A)') '# pdbfile      = output_{}.pdb   # PDB files'

      if (toks(i) == 'pdb_tgt') &
        write(MsgOut,'(A)') 'pdb_tgtfile    = output.pdb      # PDB file (Morphing)'

      if (toks(i) == 'pdb_ave') &
        write(MsgOut,'(A)') 'pdb_avefile    = output_ave.pdb  # PDB file (Averaged coordinates of analysis atoms)'

      if (toks(i) == 'pdb_aft') &
        write(MsgOut,'(A)') 'pdb_aftfile    = output_aft.pdb  # PDB file (Averaged coordinates of fitting atoms)'

      if (toks(i) == 'pmf') &
        write(MsgOut,'(A)') 'pmffile        = output.pmf      # potential of mean force file'

      if (toks(i) == 'pml') &
        write(MsgOut,'(A)') 'pmlfile        = output.pml      # PyMol script file'

      if (toks(i) == 'prj') &
        write(MsgOut,'(A)') 'prjfile        = output.prj      # PRJ file'

      if (toks(i) == 'prob') &
        write(MsgOut,'(A)') 'probfile       = output.prob     # unbiased density of states file'

      if (toks(i) == 'qmmm_crd') &
        write(MsgOut,'(A)') 'qmmm_crdfile   = snapshot{}.crd  # CHARMM CARD file for QMMM calc.'

      if (toks(i) == 'qmmm_psf') &
        write(MsgOut,'(A)') 'qmmm_psffile   = snapshot{}.psf  # CHARMM PSF file for QMMM calc.'

      if (toks(i) == 'qmmm_pdb') &
        write(MsgOut,'(A)') 'qmmm_pdbfile   = snapshot{}.pdb  # PDB file for reference of analysis'

      if (toks(i) == 'qnt') &
        write(MsgOut,'(A)') 'qntfile        = output.qnt      # fraction of native contact file'

      if (toks(i) == 'rdf') &
        write(MsgOut,'(A)') 'rdffile        = output.rdf      # radial distribution function file'

      if (toks(i) == 'rms') &
        write(MsgOut,'(A)') 'rmsfile        = output.rms      # RMSD file'

      if (toks(i) == 'rg') &
        write(MsgOut,'(A)') 'rgfile         = output.rg       # RG file'

      if (toks(i) == 'rst') &
        write(MsgOut,'(A)') 'rstfile        = output.rst      # restart file'

      if (toks(i) == 'RST') &
        write(MsgOut,'(A)') '# rstfile        = output{}.rst    # restart file'

      if (toks(i) == 'txt') &
        write(MsgOut,'(A)') 'txtfile        = output.txt      # text file'

      if (toks(i) == 'hblist') &
        write(MsgOut,'(A)') '# hb_listfile  = output.hb_list    # H-bond list file'

      if (toks(i) == 'HBLIST') &
        write(MsgOut,'(A)') '# hb_listfile  = output_().hb_list    # parallel-IO H-bond list file'

      if (toks(i) == 'top') &
        write(MsgOut,'(A)') 'topfile        = output.top      # topology file'

      if (toks(i) == 'tor') &
        write(MsgOut,'(A)') 'torfile        = output.tor      # torsion file'

      if (toks(i) == 'TRJ') &
        write(MsgOut,'(A)') 'trjfile        = output{}.trj    # trajectory file'

      if (toks(i) == 'trj') &
        write(MsgOut,'(A)') 'trjfile        = output.trj      # trajectory file'

      if (toks(i) == 'STRJ') &
        write(MsgOut,'(A)') '# trjfile        = output{}.pdb    # split PDB trajectory file'

      if (toks(i) == 'trr') &
        write(MsgOut,'(A)') 'trrfile        = output.trr      # TRROT file'

      if (toks(i) == 'val') &
        write(MsgOut,'(A)') 'valfile        = output.val      # VAL file'

      if (toks(i) == 'vcv') &
        write(MsgOut,'(A)') 'vcvfile        = output.vcv      # Variance-Covarience Matrix file'

      if (toks(i) == 'vec') &
        write(MsgOut,'(A)') 'vecfile        = output.vec      # VEC file'

      if (toks(i) == 'vel') &
        write(MsgOut,'(A)') 'velfile        = output.vel      # velocity restart file'

      if (toks(i) == 'voronoi') &
        write(MsgOut,'(A)') 'voronoifile    = output.voronoi  # vorinoi file'

      if (toks(i) == 'VORONOI') &
        write(MsgOut,'(A)') '# voronoifile    = output{}.voronoi# vorinoi file'

      if (toks(i) == 'vmd') &
        write(MsgOut,'(A)') 'vmdfile        = output.vmd      # VMD visualization state file'

      if (toks(i) == 'weight') &
        write(MsgOut,'(A)') 'weightfile     = output.weight   # weight file'

      if (toks(i) == 'WEIGHT') &
        write(MsgOut,'(A)') '# weightfile     = output{}.weight # weight file'

      if (toks(i) == 'xsc') &
        write(MsgOut,'(A)') 'xscfile        = output.xsc      # extended system file'

      if (toks(i) == 'morph') &
        write(MsgOut,'(A)') 'morphfile      = output.mor      # Morphing file'

    end do

!TM(181230)
!    if (backup_allowed()) then
!      yes_no = "yes"
!    else
!      yes_no = "no"
!    end if
!
!    write(MsgOut,'(A)') backup_control_file_option//'   = '&
!      //yes_no//'             # allow backup of existing output files'
!TM(181230)

    write(MsgOut,'(A)') ' '

    deallocate(toks)

    return

  end subroutine show_ctrl_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    read_ctrl_output
  !> @brief        read control parameters in OUTPUT section
  !! @authors      NT
  !! @param[in]    handle   : unit number of control file
  !! @param[inout] out_info : OUTPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine read_ctrl_output(handle, out_info)

    ! parameters
    character(*),            parameter :: Section = 'Output'

    ! formal argments
    integer,                 intent(in)    :: handle
    type(s_out_info),        intent(inout) :: out_info

    ! local variables
    !
    ! temporary string
    character(MaxLine)                     :: string


    ! read parameters from control file
    !
    call begin_ctrlfile_section(handle, Section)

    call read_ctrlfile_string(handle, Section, 'Ambcrdfile',  &
         out_info%ambcrdfile)
    call read_ctrlfile_string(handle, Section, 'Angfile',     &
         out_info%angfile)
    call read_ctrlfile_string(handle, Section, 'Cntfile',     &
         out_info%cntfile)
    call read_ctrlfile_string(handle, Section, 'Comangfile',  &
         out_info%comangfile)
    call read_ctrlfile_string(handle, Section, 'Comdisfile',  &
         out_info%comdisfile)
    call read_ctrlfile_string(handle, Section, 'Comtorfile',  &
         out_info%comtorfile)
    call read_ctrlfile_string(handle, Section, 'Coorfile',    &
         out_info%coorfile)
    call read_ctrlfile_string(handle, Section, 'Crdfile',     &
         out_info%crdfile)
    call read_ctrlfile_string(handle, Section, 'Crsfile',     &
         out_info%crsfile)
    call read_ctrlfile_string(handle, Section, 'Disfile',     &
         out_info%disfile)
    call read_ctrlfile_string(handle, Section, 'Enefile',     &
         out_info%enefile)
    call read_ctrlfile_string(handle, Section, 'Excfile',     &
         out_info%excfile)
    call read_ctrlfile_string(handle, Section, 'Fenefile',    &
         out_info%fenefile)
    call read_ctrlfile_string(handle, Section, 'Gprfile',     &
         out_info%gprfile)
    call read_ctrlfile_string(handle, Section, 'HB_listfile', &
         out_info%hb_listfile)
    call read_ctrlfile_string(handle, Section, 'Mapfile',     &
         out_info%mapfile)
    call read_ctrlfile_string(handle, Section, 'Msdfile',     &
         out_info%msdfile)
    call read_ctrlfile_string(handle, Section, 'Grotopfile',  &
         out_info%grotopfile)
    call read_ctrlfile_string(handle, Section, 'Grocrdfile',  &
         out_info%grocrdfile)
    call read_ctrlfile_string(handle, Section, 'Grocrd_tgtfile',  &
         out_info%grocrd_tgtfile)
    call read_ctrlfile_string(handle, Section, 'Indexfile',   &
         out_info%indexfile)
    call read_ctrlfile_string(handle, Section, 'Logfile',     &
         out_info%logfile)
    call read_ctrlfile_string(handle, Section, 'Outfile',     &
         out_info%outfile)
    call read_ctrlfile_string(handle, Section, 'Parfile',     &
         out_info%parfile)
    call read_ctrlfile_string(handle, Section, 'Pathcvfile',  &
         out_info%pathcvfile)
    call read_ctrlfile_string(handle, Section, 'Pcafile',     &
         out_info%pcafile)
    call read_ctrlfile_string(handle, Section, 'Pdbfile',     &
         out_info%pdbfile)
    call read_ctrlfile_string(handle, Section, 'Pdb_avefile', &
         out_info%pdb_avefile)
    call read_ctrlfile_string(handle, Section, 'Pdb_aftfile', &
         out_info%pdb_aftfile)
    call read_ctrlfile_string(handle, Section, 'Pmffile',     &
         out_info%pmffile)
    call read_ctrlfile_string(handle, Section, 'Pmlfile',     &
         out_info%pmlfile)
    call read_ctrlfile_string(handle, Section, 'Prjfile',     &
         out_info%prjfile)
    call read_ctrlfile_string(handle, Section, 'Probfile',    &
         out_info%probfile)
    call read_ctrlfile_string(handle, Section, 'qmmm_crdfile',&
         out_info%qmmm_crdfile)
    call read_ctrlfile_string(handle, Section, 'qmmm_psffile',&
         out_info%qmmm_psffile)
    call read_ctrlfile_string(handle, Section, 'qmmm_pdbfile',&
         out_info%qmmm_pdbfile)
    call read_ctrlfile_string(handle, Section, 'Qntfile',     &
         out_info%qntfile)
    call read_ctrlfile_string(handle, Section, 'Rdffile',     &
         out_info%rdffile)
    call read_ctrlfile_string(handle, Section, 'Rmsfile',     &
         out_info%rmsfile)
    call read_ctrlfile_string(handle, Section, 'Rgfile',      &
         out_info%rgfile)
    call read_ctrlfile_string(handle, Section, 'Rstfile',     &
         out_info%rstfile)
    call read_ctrlfile_string(handle, Section, 'Txtfile',     &
         out_info%txtfile)
    call read_ctrlfile_string(handle, Section, 'Topfile',     &
         out_info%topfile)
    call read_ctrlfile_string(handle, Section, 'Torfile',     &
         out_info%torfile)
    call read_ctrlfile_string(handle, Section, 'Trjfile',     &
         out_info%trjfile)
    call read_ctrlfile_string(handle, Section, 'Trrfile',     &
         out_info%trrfile)
    call read_ctrlfile_string(handle, Section, 'Valfile',     &
         out_info%valfile)
    call read_ctrlfile_string(handle, Section, 'Vcvfile',     &
         out_info%vcvfile)
    call read_ctrlfile_string(handle, Section, 'Vecfile',     &
         out_info%vecfile)
    call read_ctrlfile_string(handle, Section, 'Velfile',     &
         out_info%velfile)
    call read_ctrlfile_string(handle, Section, 'voronoifile', &
         out_info%voronoifile)
    call read_ctrlfile_string(handle, Section, 'Vmdfile',     &
         out_info%vmdfile)
    call read_ctrlfile_string(handle, Section, 'weightfile',  &
         out_info%weightfile)
    call read_ctrlfile_string(handle, Section, 'Xscfile',     &
         out_info%xscfile)
    call read_ctrlfile_string(handle, Section, 'Tblfile',     &
         out_info%tblfile)
    call read_ctrlfile_string(handle, Section, 'Morphfile',   &
         out_info%morphfile)

!TM(181230)
!    string = ''
!    call read_ctrlfile_string(handle, Section, backup_control_file_option,&
!         string)
!    out_info%allow_backup = parse_tristate_string(string)
!TM(181230)

    call end_ctrlfile_section(handle)

    ! write parameters to MsgOut
    !

    if (main_rank) then
      write(MsgOut,'(A)') 'Read_Ctrl_Output> Output Files'
      if (len_trim(out_info%ambcrdfile) > 0) &
        write(MsgOut,'(A20,A)') '  ambcrdfile      = ', trim(out_info%ambcrdfile)
      if (len_trim(out_info%angfile) > 0) &
        write(MsgOut,'(A20,A)') '  angfile         = ', trim(out_info%angfile)
      if (len_trim(out_info%cntfile) > 0) &
        write(MsgOut,'(A20,A)') '  cntfile         = ', trim(out_info%cntfile)
      if (len_trim(out_info%comangfile) > 0) &
        write(MsgOut,'(A20,A)') '  comangfile      = ', trim(out_info%comangfile)
      if (len_trim(out_info%comdisfile) > 0) &
        write(MsgOut,'(A20,A)') '  comdisfile      = ', trim(out_info%comdisfile)
      if (len_trim(out_info%comtorfile) > 0) &
        write(MsgOut,'(A20,A)') '  comtorfile      = ', trim(out_info%comtorfile)
      if (len_trim(out_info%coorfile) > 0) &
        write(MsgOut,'(A20,A)') '  coorfile        = ', trim(out_info%coorfile)
      if (len_trim(out_info%crdfile) > 0) &
        write(MsgOut,'(A20,A)') '  crdfile         = ', trim(out_info%crdfile)
      if (len_trim(out_info%crsfile) > 0) &
        write(MsgOut,'(A20,A)') '  crsfile         = ', trim(out_info%crsfile)
      if (len_trim(out_info%disfile) > 0) &
        write(MsgOut,'(A20,A)') '  disfile         = ', trim(out_info%disfile)
      if (len_trim(out_info%enefile) > 0) &
        write(MsgOut,'(A20,A)') '  enefile         = ', trim(out_info%enefile)
      if (len_trim(out_info%excfile) > 0) &
        write(MsgOut,'(A20,A)') '  excfile         = ', trim(out_info%excfile)
      if (len_trim(out_info%fenefile) > 0) &
        write(MsgOut,'(A20,A)') '  fenefile        = ', trim(out_info%fenefile)
      if (len_trim(out_info%gprfile) > 0) &
        write(MsgOut,'(A20,A)') '  gprfile         = ', trim(out_info%gprfile)
      if (len_trim(out_info%hb_listfile) > 0) &
        write(MsgOut,'(A20,A)') '  hb_listfile     = ', trim(out_info%hb_listfile)
      if (len_trim(out_info%mapfile) > 0) &
        write(MsgOut,'(A20,A)') '  mapfile         = ', trim(out_info%mapfile)
      if (len_trim(out_info%msdfile) > 0) &
        write(MsgOut,'(A20,A)') '  msdfile         = ', trim(out_info%msdfile)
      if (len_trim(out_info%grotopfile) > 0) &
        write(MsgOut,'(A20,A)') '  grotopfile      = ', trim(out_info%grotopfile)
      if (len_trim(out_info%grocrdfile) > 0) &
        write(MsgOut,'(A20,A)') '  grocrdfile      = ', trim(out_info%grocrdfile)
      if (len_trim(out_info%grocrd_tgtfile) > 0) &
        write(MsgOut,'(A20,A)') '  grocrd_tgtfile  = ', &
                                                   trim(out_info%grocrd_tgtfile)
      if (len_trim(out_info%indexfile) > 0) &
        write(MsgOut,'(A20,A)') '  indexfile       = ', trim(out_info%indexfile)
      if (len_trim(out_info%logfile) > 0) &
        write(MsgOut,'(A20,A)') '  logfile         = ', trim(out_info%logfile)
      if (len_trim(out_info%outfile) > 0) &
        write(MsgOut,'(A20,A)') '  outfile         = ', trim(out_info%outfile)
      if (len_trim(out_info%parfile) > 0) &
        write(MsgOut,'(A20,A)') '  parfile         = ', trim(out_info%parfile)
      if (len_trim(out_info%pathcvfile) > 0) &
        write(MsgOut,'(A20,A)') '  pathcvfile      = ', trim(out_info%pathcvfile)
      if (len_trim(out_info%pcafile) > 0) &
        write(MsgOut,'(A20,A)') '  pcafile         = ', trim(out_info%pcafile)
      if (len_trim(out_info%pdbfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdbfile         = ', trim(out_info%pdbfile)
      if (len_trim(out_info%pdb_tgtfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_tgtfile     = ', trim(out_info%pdb_tgtfile)
      if (len_trim(out_info%pdb_avefile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_avefile     = ', trim(out_info%pdb_avefile)
      if (len_trim(out_info%pdb_aftfile) > 0) &
        write(MsgOut,'(A20,A)') '  pdb_aftfile     = ', trim(out_info%pdb_aftfile)
      if (len_trim(out_info%pmffile) > 0) &
        write(MsgOut,'(A20,A)') '  pmffile         = ', trim(out_info%pmffile)
      if (len_trim(out_info%pmlfile) > 0) &
        write(MsgOut,'(A20,A)') '  pmlfile         = ', trim(out_info%pmlfile)
      if (len_trim(out_info%prjfile) > 0) &
        write(MsgOut,'(A20,A)') '  prjfile         = ', trim(out_info%prjfile)
      if (len_trim(out_info%probfile) > 0) &
        write(MsgOut,'(A20,A)') '  probfile        = ', trim(out_info%probfile)
      if (len_trim(out_info%qmmm_crdfile) > 0) &
        write(MsgOut,'(A20,A)') '  qmmm_crdfile    = ', trim(out_info%qmmm_crdfile)
      if (len_trim(out_info%qmmm_psffile) > 0) &
        write(MsgOut,'(A20,A)') '  qmmm_psffile    = ', trim(out_info%qmmm_psffile)
      if (len_trim(out_info%qmmm_pdbfile) > 0) &
        write(MsgOut,'(A20,A)') '  qmmm_pdbfile    = ', trim(out_info%qmmm_pdbfile)
      if (len_trim(out_info%qntfile) > 0) &
        write(MsgOut,'(A20,A)') '  qntfile         = ', trim(out_info%qntfile)
      if (len_trim(out_info%rdffile) > 0) &
        write(MsgOut,'(A20,A)') '  rdffile         = ', trim(out_info%rdffile)
      if (len_trim(out_info%rmsfile) > 0) &
        write(MsgOut,'(A20,A)') '  rmsfile         = ', trim(out_info%rmsfile)
      if (len_trim(out_info%rgfile) > 0) &
        write(MsgOut,'(A20,A)') '  rgfile          = ', trim(out_info%rgfile)
      if (len_trim(out_info%rstfile) > 0) &
        write(MsgOut,'(A20,A)') '  rstfile         = ', trim(out_info%rstfile)
      if (len_trim(out_info%txtfile) > 0) &
        write(MsgOut,'(A20,A)') '  txtfile         = ', trim(out_info%txtfile)
      if (len_trim(out_info%topfile) > 0) &
        write(MsgOut,'(A20,A)') '  topfile         = ', trim(out_info%topfile)
      if (len_trim(out_info%torfile) > 0) &
        write(MsgOut,'(A20,A)') '  torfile         = ', trim(out_info%torfile)
      if (len_trim(out_info%trjfile) > 0) &
        write(MsgOut,'(A20,A)') '  trjfile         = ', trim(out_info%trjfile)
      if (len_trim(out_info%trrfile) > 0) &
        write(MsgOut,'(A20,A)') '  trrfile         = ', trim(out_info%trrfile)
      if (len_trim(out_info%valfile) > 0) &
        write(MsgOut,'(A20,A)') '  valfile         = ', trim(out_info%valfile)
      if (len_trim(out_info%vcvfile) > 0) &
        write(MsgOut,'(A20,A)') '  vcvfile         = ', trim(out_info%vcvfile)
      if (len_trim(out_info%vecfile) > 0) &
        write(MsgOut,'(A20,A)') '  vecfile         = ', trim(out_info%vecfile)
      if (len_trim(out_info%velfile) > 0) &
        write(MsgOut,'(A20,A)') '  velfile         = ', trim(out_info%velfile)
      if (len_trim(out_info%vmdfile) > 0) &
        write(MsgOut,'(A20,A)') '  vmdfile         = ', trim(out_info%vmdfile)
      if (len_trim(out_info%weightfile) > 0) &
        write(MsgOut,'(A20,A)') '  weightfile      = ', trim(out_info%weightfile)
      if (len_trim(out_info%xscfile) > 0) &
        write(MsgOut,'(A20,A)') '  xscfile         = ', trim(out_info%xscfile)
      if (len_trim(out_info%tblfile) > 0) &
        write(MsgOut,'(A20,A)') '  tblfile         = ', trim(out_info%tblfile)
      if (len_trim(out_info%morphfile) > 0) &
        write(MsgOut,'(A20,A)') '  morphfile       = ', trim(out_info%morphfile)

!TM(181230)
!      select case(out_info%allow_backup)
!      case(tristate_TRUE)
!        string = "yes"
!      case(tristate_FALSE)
!        string = "no"
!      case(tristate_NOT_SET)
!        string = "not set"
!      end select
!      write(MsgOut,'(A20,A)') '  '//backup_control_file_option//&
!        '    = ', trim(string)
!TM(181230)

      write(MsgOut,'(A)') ' '

    end if
 
    return

  end subroutine read_ctrl_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    parse_ctrl_output
  !> @brief        parse control parameters in OUTPUT section from commandline
  !! @authors      MK
  !! @param[in]    section  : output section data
  !! @param[inout] out_info : OUTPUT section control parameters information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine parse_ctrl_output( section, out_info )

    type(genesis_opt_section), intent(in)    :: section
    type(s_out_info),          intent(inout) :: out_info

    integer :: i

    do i = 1, section%count
      select case( section%items(i)%name )
        case ("angfile")
          out_info%angfile = section%items(i)%value
        case ("cntfile")
          out_info%cntfile = section%items(i)%value
        case ("comangfile")
          out_info%comangfile = section%items(i)%value
        case ("comdisfile")
          out_info%comdisfile = section%items(i)%value
        case ("comtorfile")
          out_info%comtorfile = section%items(i)%value
        case ("coorfile")
          out_info%coorfile = section%items(i)%value
        case ("crsfile")
          out_info%crsfile = section%items(i)%value
        case ("disfile")
          out_info%disfile = section%items(i)%value
        case ("enefile")
          out_info%enefile = section%items(i)%value
        case ("fenefile")
          out_info%fenefile = section%items(i)%value
        case ("gprfile")
          out_info%gprfile = section%items(i)%value
        case ("grotopfile")
          out_info%grotopfile = section%items(i)%value
        case ("grocrdfile")
          out_info%grocrdfile = section%items(i)%value
        case ("grocrd_tgtfile")
          out_info%grocrd_tgtfile = section%items(i)%value
        case ("pathcvfile")
          out_info%pathcvfile = section%items(i)%value
        case ("pcafile")
          out_info%pcafile = section%items(i)%value
        case ("pdbfile")
          out_info%pdbfile = section%items(i)%value
        case ("pdb_avefile")
          out_info%pdb_avefile = section%items(i)%value
        case ("pdb_aftfile")
          out_info%pdb_aftfile = section%items(i)%value
        case ("pmffile")
          out_info%pmffile = section%items(i)%value
        case ("prjfile")
          out_info%prjfile = section%items(i)%value
        case ("probfile")
          out_info%probfile = section%items(i)%value
        case ("qntfile")
          out_info%qntfile = section%items(i)%value
        case ("rmsfile")
          out_info%rmsfile = section%items(i)%value
        case ("rstfile")
          out_info%rstfile = section%items(i)%value
        case ("torfile")
          out_info%torfile = section%items(i)%value
        case ("trjfile")
          out_info%trjfile = section%items(i)%value
        case ("trrfile")
          out_info%trrfile = section%items(i)%value
        case ("valfile")
          out_info%valfile = section%items(i)%value
        case ("vcvfile")
          out_info%vcvfile = section%items(i)%value
        case ("vecfile")
          out_info%vecfile = section%items(i)%value
        case ("velfile")
          out_info%velfile = section%items(i)%value
        case ("voronoifile")
          out_info%voronoifile = section%items(i)%value
        case ("weightfile")
          out_info%weightfile = section%items(i)%value
        case ("xscfile")
          out_info%xscfile = section%items(i)%value
        case ("morphfile")
          out_info%morphfile = section%items(i)%value
        case( "indexfile")
          out_info%indexfile = section%items(i)%value
        case default
          write(MsgOut,'("Error(OUTPUT)> unknown param name ", a, " found.")') &
              trim(section%items(i)%name)
          call error_msg('Error(parse_ctrl_output)> quit program')
      end select

      write(MsgOut,'("INFO> PARAMETER:", a, " of SECTION:", &
                          & a, " is replaced to:", a )') &
          trim(section%items(i)%name), trim(section%name), &
          trim(section%items(i)%value)

    end do

    return

  end subroutine parse_ctrl_output

  !======1=========2=========3=========4=========5=========6=========7=========8
  !
  !  Subroutine    setup_output
  !> @brief        setup output files
  !! @authors      NT
  !! @param[in]    out_info : OUTPUT section control parameters information
  !! @param[out]   output   : output information
  !
  !======1=========2=========3=========4=========5=========6=========7=========8

  subroutine setup_output(out_info, output)

    ! formal argments
    type(s_out_info),        intent(in)    :: out_info
    type(s_output),          intent(inout) :: output


    output%ambcrdfile   = out_info%ambcrdfile
    output%angfile      = out_info%angfile
    output%cntfile      = out_info%cntfile
    output%comangfile   = out_info%comangfile
    output%comdisfile   = out_info%comdisfile
    output%comtorfile   = out_info%comtorfile
    output%coorfile     = out_info%coorfile
    output%crdfile      = out_info%crdfile
    output%crsfile      = out_info%crsfile
    output%disfile      = out_info%disfile
    output%enefile      = out_info%enefile
    output%excfile      = out_info%excfile
    output%fenefile     = out_info%fenefile
    output%gprfile      = out_info%gprfile
    output%grotopfile   = out_info%grotopfile
    output%grocrdfile   = out_info%grocrdfile
    output%grocrd_tgtfile  = out_info%grocrd_tgtfile
    output%hb_listfile  = out_info%hb_listfile
    output%indexfile    = out_info%indexfile
    output%logfile      = out_info%logfile
    output%mapfile      = out_info%mapfile
    output%msdfile      = out_info%msdfile
    output%outfile      = out_info%outfile
    output%parfile      = out_info%parfile
    output%pathcvfile   = out_info%pathcvfile
    output%pcafile      = out_info%pcafile
    output%pdbfile      = out_info%pdbfile
    output%pdb_avefile  = out_info%pdb_avefile
    output%pdb_aftfile  = out_info%pdb_aftfile
    output%pdb_tgtfile  = out_info%pdb_tgtfile
    output%pmffile      = out_info%pmffile
    output%pmlfile      = out_info%pmlfile
    output%prjfile      = out_info%prjfile
    output%probfile     = out_info%probfile
    output%qmmm_crdfile = out_info%qmmm_crdfile
    output%qmmm_psffile = out_info%qmmm_psffile
    output%qmmm_pdbfile = out_info%qmmm_pdbfile
    output%qntfile      = out_info%qntfile
    output%rdffile      = out_info%rdffile
    output%rmsfile      = out_info%rmsfile
    output%rgfile       = out_info%rgfile
    output%rstfile      = out_info%rstfile
    output%txtfile      = out_info%txtfile
    output%topfile      = out_info%topfile
    output%torfile      = out_info%torfile
    output%trjfile      = out_info%trjfile
    output%trrfile      = out_info%trrfile
    output%valfile      = out_info%valfile
    output%vcvfile      = out_info%vcvfile
    output%vecfile      = out_info%vecfile
    output%velfile      = out_info%velfile
    output%vmdfile      = out_info%vmdfile
    output%weightfile   = out_info%weightfile
    output%xscfile      = out_info%xscfile
    output%tblfile      = out_info%tblfile
    output%morphfile    = out_info%morphfile

!TM(181230)
!    call set_backup_allowed_by_control(out_info%allow_backup)
!TM(181230)

    return

  end subroutine setup_output

end module output_mod
