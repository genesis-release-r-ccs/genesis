!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   input_str_mod
!> @brief   structure of input information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

module input_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_input
    character(MaxFilename) :: ambcrdfile
    character(MaxFilename) :: ambreffile
    character(MaxFilename) :: atpfile
    character(MaxFilename) :: coorfile
    character(MaxFilename) :: cvfile
    character(MaxFilename) :: targetfile
    character(MaxFilename) :: dcdfile
    character(MaxFilename) :: dcdvelfile
    character(MaxFilename) :: excfile
    character(MaxFilename) :: gprfile
    character(MaxFilename) :: grocrdfile
    character(MaxFilename) :: groreffile
    character(MaxFilename) :: grotopfile
    character(MaxFilename) :: mtfile
    character(MaxFilename) :: indexfile
    character(MaxFilename) :: logfile
    character(MaxFilename) :: enefile
    character(MaxFilename) :: msdfile
    character(MaxFilename) :: pathfile
    character(MaxFilename) :: pathcvfile
    character(MaxFilename) :: pcafile
    character(MaxFilename) :: pdbfile
    character(MaxFilename) :: pdb_tgtfile
    character(MaxFilename) :: pdb_avefile
    character(MaxFilename) :: pdb_aftfile
    character(MaxFilename) :: pdb_sphfile
    character(MaxFilename) :: pdb_wbxfile
    character(MaxFilename) :: prmtopfile
    character(MaxFilename) :: psffile
    character(MaxFilename) :: radfile
    character(MaxFilename) :: refenefile
    character(MaxFilename) :: reffile
    character(MaxFilename) :: fitfile
    character(MaxFilename) :: remfile
    character(MaxFilename) :: rstfile
    character(MaxFilename) :: rtpfile
    character(MaxFilename) :: topfile
    character(MaxFilename) :: valfile
    character(MaxFilename) :: vecfile
    character(MaxFilename) :: velfile
    character(MaxFilename) :: weightfile
    character(MaxFilename) :: xscfile
    character(MaxFilename) :: distfile
  end type s_input

end module input_str_mod
