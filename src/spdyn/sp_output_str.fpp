!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   sp_output_str_md
!> @brief   structure of output information
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module sp_output_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_output
    integer                :: logunit = MsgOut
    integer                :: dcdunit
    integer                :: dcdvelunit
    integer                :: eneunit
    integer                :: rstunit
    integer                :: pdbunit
    integer                :: remunit
    integer                :: rpathunit
    integer                :: rpathlogunit
    integer                :: mfrcunit
    integer                :: gamdunit
    integer                :: fepunit

    character(MaxFilename) :: logfile
    character(MaxFilename) :: dcdfile
    character(MaxFilename) :: dcdvelfile
    character(MaxFilename) :: enefile
    character(MaxFilename) :: rstfile
    character(MaxFilename) :: pdbfile
    character(MaxFilename) :: remfile
    character(MaxFilename) :: rpathfile
    character(MaxFilename) :: rpathlogfile
    character(MaxFilename) :: mfrcfile
    character(MaxFilename) :: gamdfile
    character(MaxFilename) :: fepfile

    logical                :: out_dcd_header    = .true.
    logical                :: out_dcdvel_header = .true.

    logical                :: logout     = .false.
    logical                :: dcdout     = .false.
    logical                :: dcdvelout  = .false.
    logical                :: eneout     = .false.
    logical                :: rstout     = .false.
    logical                :: pdbout     = .false.
    logical                :: remout     = .false.
    logical                :: rpathout   = .false.
    logical                :: rpathlogout= .false.
    logical                :: mfrcout    = .false.
    logical                :: gamdout    = .false.
    logical                :: fepout     = .false.
                    
    logical                :: replica    = .false.
    logical                :: rpath      = .false.

    logical                :: dcd_para
    logical                :: dcdvel_para
    logical                :: rst_para

    logical                :: verbose

  end type s_output

end module sp_output_str_mod
