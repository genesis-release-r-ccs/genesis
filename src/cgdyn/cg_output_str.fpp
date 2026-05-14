!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   cg_output_str_md
!> @brief   Defines the structure of output information.
!! @authors Takaharu Mori (TM)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

module cg_output_str_mod

  use string_mod
  use messages_mod

  implicit none
  private

  ! structures
  type, public :: s_output
    integer                :: logunit = MsgOut
    integer                :: dcdunit
    integer                :: dcdvelunit
    integer                :: rstunit
    integer                :: pdbunit
    integer                :: remunit
    integer                :: mfrcunit

    character(MaxFilename) :: logfile
    character(MaxFilename) :: dcdfile
    character(MaxFilename) :: dcdvelfile
    character(MaxFilename) :: rstfile
    character(MaxFilename) :: pdbfile
    character(MaxFilename) :: remfile
    character(MaxFilename) :: mfrcfile

    logical                :: out_dcd_header    = .true.
    logical                :: out_dcdvel_header = .true.

    logical                :: logout     = .false.
    logical                :: dcdout     = .false.
    logical                :: dcdvelout  = .false.
    logical                :: rstout     = .false.
    logical                :: pdbout     = .false.
    logical                :: remout     = .false.
    logical                :: mfrcout    = .false.
                    
    logical                :: replica    = .false.

    logical                :: dcd_para
    logical                :: dcdvel_para
    logical                :: rst_para

    logical                :: verbose

  end type s_output

end module cg_output_str_mod
