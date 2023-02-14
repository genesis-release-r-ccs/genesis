#!/bin/bash

# This script cleans up output and log files created 
# by regression tests, and restores the status before
# the tests.
#
# Usage:
#    $ ./cleanup.sh

rm test_atdyn/*/*/error*     >& /dev/null
rm test_atdyn/*/*/test*      >& /dev/null
rm -r test_atdyn/*/*/output.?*      >& /dev/null
rm test_atdyn/br/*/qmmm.0/*.inp            >& /dev/null
rm test_atdyn/br/*/qmmm.0/*qm.xyz          >& /dev/null
rm test_atdyn/br/*/qmmm.0/*pc.xyz          >& /dev/null
rm test_atdyn/br/*/qmmm.0/mm_charges.dat   >& /dev/null

rm test_spdyn/*/*/error*     >& /dev/null
rm test_spdyn/*/*/test*      >& /dev/null
rm -r test_spdyn/*/*/output.?*      >& /dev/null

rm test_remd_atdyn/*/error*  >& /dev/null
rm test_remd_atdyn/*/log     >& /dev/null
rm test_remd_atdyn/*/test?   >& /dev/null
rm -r test_remd_atdyn/*/output.?*      >& /dev/null
rm test_remd_spdyn/*/error*  >& /dev/null
rm test_remd_spdyn/*/log     >& /dev/null
rm test_remd_spdyn/*/test?   >& /dev/null
rm -r test_remd_spdyn/*/output.?*      >& /dev/null

rm test_rpath_atdyn/*/error* >& /dev/null
rm test_rpath_atdyn/*/log    >& /dev/null
rm test_rpath_atdyn/*/test*  >& /dev/null
rm test_rpath_atdyn/tim*/dcd*          >& /dev/null
rm test_rpath_atdyn/tim*/rst*          >& /dev/null
rm test_rpath_atdyn/tim*/qmmm.*/*.inp           >& /dev/null
rm test_rpath_atdyn/tim*/qmmm.*/mm_charges.dat  >& /dev/null
rm -r test_rpath_atdyn/*/output.?*      >& /dev/null

rm test_rpath_spdyn/*/error* >& /dev/null
rm test_rpath_spdyn/*/log    >& /dev/null
rm test_rpath_spdyn/*/test?  >& /dev/null
rm -r test_rpath_spdyn/*/output.?*      >& /dev/null

rm test_vib/*/error*         >& /dev/null
rm test_vib/*/log            >& /dev/null
rm test_vib/*/vib.minfo      >& /dev/null
rm test_vib/*/test?          >& /dev/null
rm test_vib/*/qmmm.?/*.inp   >& /dev/null
rm -r test_vib/*/minfo.files >& /dev/null
rm -r test_vib/*/output.?*      >& /dev/null

rm test_gamd_atdyn/*/*/test*    >& /dev/null
rm test_gamd_atdyn/*/*/error*   >& /dev/null
rm test_gamd_atdyn/*/*/out.gamd >& /dev/null
rm -r test_gamd_atdyn/*/*/output.?*      >& /dev/null

rm test_gamd_spdyn/*/*/test*    >& /dev/null
rm test_gamd_spdyn/*/*/error*   >& /dev/null
rm test_gamd_spdyn/*/*/out.gamd >& /dev/null
rm -r test_gamd_spdyn/*/*/output.?*      >& /dev/null

rm test_fep/*/*/test*       >& /dev/null
rm test_fep/*/*/error*      >& /dev/null
rm test_fep/*/*/out*.fepout >& /dev/null
rm -r test_fep/*/*/output.?*      >& /dev/null

rm test_parallel_IO/*/*/error*     >& /dev/null
rm test_parallel_IO/*/*/test*      >& /dev/null
rm test_parallel_IO/*/*/log*      >& /dev/null
rm test_parallel_IO/*/*/ref*      >& /dev/null
rm test_parallel_IO/*/*/*.rst      >& /dev/null
rm -r test_parallel_IO/*/*/cache      >& /dev/null

cd test_analysis
./cleanup.sh  >& /dev/null
cd ..
cd test_spana
./cleanup.sh  >& /dev/null
cd ..
