#!/usr/bin/env python
# coding: utf-8

import commands
import os
import os.path
import sys
import copy
import random
from genesis import *

############### MAIN ##########################################################

## configuration

# default tests
test_system_common = [
'dna',
'dppc',
'jac',
]

test_list_common = [
'LEAP_CUTOFF'              ,
'VVER_CUTOFF'              ,
'LEAP_CUTOFF_FSW'          ,
'LEAP_CUTOFF_REST-ANGLE'   ,
'LEAP_CUTOFF_REST-DIHED'   ,
'LEAP_CUTOFF_REST-DIST'    ,
'LEAP_CUTOFF_REST-POSI'    ,
'LEAP_PME'                 ,
'LEAP_SHAKE_CUTOFF'        ,
'VVER_RATTLE_CUTOFF'       ,
'LEAP_SHAKE_PME',
'LEAP_SHAKE_FSW_TIP3',
'LEAP_SHAKE_PME_TABLE-1',
'LEAP_SHAKE_PME_TABLE-1_FSW',
'LEAP_SHAKE_PME_TABLE-3',
'LEAP_SHAKE_PME_TABLE-3_FSW',
'LEAP_SHAKE_PME_TABLE-3_FSW_TIP3',
'LEAP_SHAKE_PME_TABLE-3_TIP3',
'LEAP_SHAKE_PME_TIP3',
]

test_list_spdyn = [
'LEAP_CUTOFF_REST-LANGLE'  ,
'LEAP_CUTOFF_REST-LDIHED'  ,
'LEAP_CUTOFF_REST-LDIST'   ,
]

test_list_ensemble = [
'LEAP_NPT_SHAKE_BERENDSEN' ,
'LEAP_NVT_SHAKE_BERENDSEN' ,
'VVER_NPT_RATTLE_LANGEVIN'  ,
'LEAP_NPT_SHAKE_LANGEVIN'  ,
'LEAP_NVT_SHAKE_LANGEVIN'  ,
'VVER_NVT_RATTLE_LANGEVIN'  ,
]

test_list_atdyn = [
'LEAP_CUTOFF_REST-RMSD'    , 
]

# global variables
tolerance_strict = 0.0001 # energy difference in kcal/mol
#os.environ["OMP_NUM_THREADS"] = "1"
ipassed = 0
ifailed = 0
iaborted = 0
itried = 0

## parse command line
if len(sys.argv) == 1:
    genesis_command = 'mpirun -np 8 atdyn'
elif len(sys.argv) == 2:
    genesis_command = sys.argv[1]
else:
    genesis_command = sys.argv[1]
    test_list = sys.argv[2:]

if (genesis_command == "-h") or (genesis_command == "--help"):
    print """    usage:
    $ ./test.py [genesis command] [directory names]
    
    examples:
    # run tests using the default command and directories
    $ ./test.py

    # run atdyn
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/atdyn"

    # run spdyn
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/spdyn"

    # run parallel IO
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/spdyn" parallel_io
    """
    sys.exit(3)

genesis_command_last = genesis_command.split()[-1]
if not os.path.exists(os.path.expanduser(genesis_command_last)):
    print "Error: %s does not exist" % genesis_command_last
    sys.exit(3)

test_list_md          = []
test_list_parallel    = [] #todo

for system_each in test_system_common:
    for test_each in test_list_common:
        dirname = 'test_common/%s/%s' % (system_each, test_each)
        if not os.path.isdir(dirname): 
            continue
        test_list_md.append(dirname)
 
pattern_parallel = re.compile('^parall')
result = None
if len(sys.argv) == 3:
    result = pattern_parallel.search(sys.argv[2])

if genesis_command.find('spdyn') == -1:
    #atdyn
    for test_each in test_list_atdyn:
        dirname = 'test_atdyn/jac/%s' % test_each
        if not os.path.isdir(dirname): 
            continue
        test_list_md.append(dirname)
    
    for test_each in test_list_ensemble:
        dirname = 'test_atdyn/dppc/%s' % test_each
        if not os.path.isdir(dirname): 
            continue
        test_list_md.append(dirname)
elif result is not None:
    #spdyn parallel IO
    for test_each in test_list_parallel:
        dirname = 'test_parallelIO/%s' % test_each
        if not os.path.isdir(dirname): 
            continue
        test_list_parallel.append(dirname)
else:
    #spdyn
    for test_each in test_list_spdyn:
        dirname = 'test_spdyn/jac/%s' % test_each
        if not os.path.isdir(dirname): 
            continue
        test_list_md.append(dirname)
    
    for test_each in test_list_ensemble:
        dirname = 'test_spdyn/dppc/%s' % test_each
        if not os.path.isdir(dirname): 
            continue
        test_list_md.append(dirname)

## common tests
if len(test_list_md):
    print "======================================================================="
    print " Regression tests for MD"
    print "======================================================================="

cwdname = os.getcwd()
for test_each in test_list_md:
    itried = itried + 1
    os.chdir(cwdname)
    dirname = test_each
    if not os.path.isdir(dirname) : 
        continue
    os.chdir(dirname)

    # run MD
    print "-----------------------------------------------------------------------"
    print "Running %s..." % (dirname + "/")

    inputname = "inp"
    testname = "test"
    commandline = '%s %s 1> %s 2> error' % (genesis_command, inputname, testname)
    print "$ %s" % commandline
    status = commands.getstatusoutput(commandline)

    if status[0] > 0:
        print
        print "Aborted..."
        print
        iaborted = iaborted + 1
        continue

    # parse the result
    #if status[0] == 0:
    ref = Genesis()
    refname = "ref"
    ref.read(refname)
    test = Genesis()
    test.read(testname)

    # check the result
    print
    print "Checking %s" % test_each
    print "Checking diff between %s and %s..." % (refname, testname)
    print
    #ref.test_diff_energies(test, tolerance_strict)
    ref.test_diff(test, tolerance_strict)
    print
    if ref.is_passed:
        ipassed = ipassed + 1
    else:
        ifailed = ifailed + 1
    # post-cleaning
    #os.remove(testname)

if (itried > 0):
    print "-----------------------------------------------------------------------"
    print "Passed  %d / %d" % (ipassed,  itried)
    print "Failed  %d / %d" % (ifailed,  itried)
    print "Aborted %d / %d" % (iaborted, itried)
    print "-----------------------------------------------------------------------"

if iaborted > 0:
    sys.exit(2)
elif ifailed > 0:
    sys.exit(1)
else:
    sys.exit(0)

