!--------1---------2---------3---------4---------5---------6---------7---------8
!
!  Module   pr_huge_molecule_mod
!> @brief   huge molecule information
!! @authors Norio Takase (NT)
!
!  (c) Copyright 2014 RIKEN. All rights reserved.
!
!--------1---------2---------3---------4---------5---------6---------7---------8

#ifdef HAVE_CONFIG_H
#include "../config.h"
#endif

#ifdef HM_DISK
#include "pr_huge_molecule_file.fpp"
#else
#include "pr_huge_molecule_mem.fpp"
#endif
