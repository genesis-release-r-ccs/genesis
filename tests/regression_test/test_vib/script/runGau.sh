#!/bin/bash

# -----------------------------------------------
# restart from given log and Fchk file
# -----------------------------------------------

QMINP=$1
QMOUT=$2
MOL=${QMINP%.*}

if [ -e ${MOL}.Fchk ]; then
   cp ${MOL}.Fchk gaussian.Fchk
   exit 0
fi

