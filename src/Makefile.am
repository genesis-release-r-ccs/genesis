#--------1---------2---------3---------4---------5---------6---------7---------8
#
#  File     : Makefile.am
#  Authors  : NT
#
#  (c) Copyright 2014 RIKEN. All rights reserved.
#
#--------1---------2---------3---------4---------5---------6---------7---------8

if ATDYN
  ATDYN_DIR=atdyn
else
  ATDYN_DIR=
endif
if ANALYSIS
  ANALYSIS_DIR=analysis
else
  ANALYSIS_DIR=
endif
if SPDYN
  SPDYN_DIR=spdyn
else
  SPDYN_DIR=
endif
SUBDIRS = lib $(ATDYN_DIR) $(SPDYN_DIR) $(ANALYSIS_DIR)

MAINTAINERCLEANFILES = Makefile.in

EXTRA_DIST = \
	GENESIS_VERSION

