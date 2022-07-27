!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   output_str_mod
!> @brief   structure of output information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module output_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_output
    character(MaxFilename) :: ambcrdfile
    character(MaxFilename) :: angfile
    character(MaxFilename) :: cntfile
    character(MaxFilename) :: comangfile
    character(MaxFilename) :: comdisfile
    character(MaxFilename) :: comtorfile
    character(MaxFilename) :: coorfile
    character(MaxFilename) :: crdfile
    character(MaxFilename) :: crsfile
    character(MaxFilename) :: disfile
    character(MaxFilename) :: enefile
    character(MaxFilename) :: excfile
    character(MaxFilename) :: fenefile
    character(MaxFilename) :: gprfile
    character(MaxFilename) :: grotopfile
    character(MaxFilename) :: grocrdfile
    character(MaxFilename) :: grocrd_tgtfile
    character(MaxFilename) :: hb_listfile
    character(MaxFilename) :: indexfile
    character(MaxFilename) :: logfile
    character(MaxFilename) :: mapfile
    character(MaxFilename) :: morphfile
    character(MaxFilename) :: msdfile
    character(MaxFilename) :: outfile
    character(MaxFilename) :: parfile
    character(MaxFilename) :: pathcvfile
    character(MaxFilename) :: pcafile
    character(MaxFilename) :: pdbfile
    character(MaxFilename) :: pdb_tgtfile
    character(MaxFilename) :: pdb_avefile
    character(MaxFilename) :: pdb_aftfile
    character(MaxFilename) :: pmffile
    character(MaxFilename) :: pmlfile
    character(MaxFilename) :: prjfile
    character(MaxFilename) :: probfile
    character(MaxFilename) :: qmmm_crdfile
    character(MaxFilename) :: qmmm_psffile
    character(MaxFilename) :: qmmm_pdbfile
    character(MaxFilename) :: qntfile
    character(MaxFilename) :: rdffile
    character(MaxFilename) :: rmsfile
    character(MaxFilename) :: rgfile
    character(MaxFilename) :: rstfile
    character(MaxFilename) :: txtfile
    character(MaxFilename) :: tblfile
    character(MaxFilename) :: topfile
    character(MaxFilename) :: torfile
    character(MaxFilename) :: trjfile
    character(MaxFilename) :: trrfile
    character(MaxFilename) :: valfile
    character(MaxFilename) :: vcvfile
    character(MaxFilename) :: vecfile
    character(MaxFilename) :: velfile
    character(MaxFilename) :: voronoifile
    character(MaxFilename) :: weightfile
    character(MaxFilename) :: xscfile
    character(MaxFilename) :: vmdfile
  end type s_output

end module output_str_mod
