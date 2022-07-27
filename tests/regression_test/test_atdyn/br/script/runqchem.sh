#!/bin/bash

# -----------------------------------------------
# restart from given log and Fchk file
# -----------------------------------------------

QMINP=$1
QMOUT=$2
MOL=${QMINP%.*}

if [ -e ${MOL}_efield.dat ]; then
   cp ${MOL}_efield.dat efield.dat
fi

