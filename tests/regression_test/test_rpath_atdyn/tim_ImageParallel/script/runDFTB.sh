#!/bin/bash

# -----------------------------------------------
# restart from given log and Fchk file
# -----------------------------------------------

QMINP=$1
QMOUT=$2
MOL=${QMINP%.*}

if [ -e ${MOL}_charges.bin ]; then
   cp ${MOL}_charges.bin  charges.bin
   cp ${MOL}_detailed.out detailed.out
   exit 0
fi

