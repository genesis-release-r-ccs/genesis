#!/usr/bin/env python
# coding: utf-8

#
# A python script for GENESIS regression tests
#
# (c) Copyright 2022 RIKEN. All rights reserved.
#

import subprocess
import os
import os.path
import sys
import copy
import random
import glob
import shutil
import re
from genesis import *


############### DEFINITION ##################################
def getdirs(path):
    test_dirs = []
    for forcefield_dir in os.listdir(path):
        forcefield_dir_path = os.path.join(path,forcefield_dir)
        if os.path.isdir(forcefield_dir_path):
            for test_dir in os.listdir(forcefield_dir_path):
                test_dir_path = os.path.join(forcefield_dir_path,test_dir)
                if os.path.isdir(test_dir_path):
                    if os.path.exists(test_dir_path + "/inp") and os.path.exists(test_dir_path + "/ref"):
                        test_dirs.append(test_dir_path)
    test_dirs.sort()
    return test_dirs

############### MAIN ########################################

###### initialization

#os.environ["OMP_NUM_THREADS"] = "1"
tolerance = 1.0e-8 # relative energy difference (diff/abs(e))
tolerance_fujitsu = 5.0e-6 # relative energy difference (diff/abs(e))
tolerance_gpu_respa = 5.0e-6 # relative energy difference (diff/abs(e))
tolerance_fujitsu_weak = 1.0e-5 # relative energy difference (diff/abs(e))
tolerance_single = 3.0e-5 # relative energy difference (diff/abs(e))
virial_ratio = 1.0e2

ipassed = 0
ifailed = 0
iaborted = 0
itried = 0
num = 0

is_atdyn = False
is_spdyn = False
is_parallelio = False
is_gpu = False
is_fugaku = False

test_dirs = []

###### parse command line

if len(sys.argv) == 1:
    genesis_command = 'mpirun -np 8 atdyn'
elif len(sys.argv) == 2:
    genesis_command = sys.argv[1]
else:
    genesis_command = sys.argv[1]
    is_number = re.compile('^[+-]?(\d*\.\d+|\d+\.?\d*)([eE][+-]?\d+|)\Z')
    if is_number.match(sys.argv[2]):
        tolerance = float(sys.argv[2])
    elif sys.argv[2] == "parallel_io":
        is_parallelio = True
    elif sys.argv[2] == "fugaku":
        mpiexec_command = genesis_command.split(' ',2)[0]
        genesis_path = genesis_command.split(' ',2)[1]
        is_fugaku = True
        if len(sys.argv) == 4:
            if is_number.match(sys.argv[3]):
                tolerance = float(sys.argv[3])
    else:
        test_dirs = sys.argv[2:]

if (genesis_command == "-h") or (genesis_command == "--help"):
    print("""    usage:
    $ ./test.py ["genesis command"] [parallel_io or tolerance_value or directories]
    
    examples:
    # run tests using the default command ("mpirun -np 8 atdyn")
    $ ./test_gamd.py

    # run atdyn tests
    $ ./test_gamd.py "mpirun -np 8 /path/to/genesis/bin/atdyn"

    # run spdyn tests
    $ ./test_gamd.py "mpirun -np 8 /path/to/genesis/bin/spdyn"

    # run with specfic tolerance value
    $ ./test_gamd.py "mpirun -np 8 /path/to/genesis/bin/atdyn" 0.1

    """)
    sys.exit(3)

genesis_command_split = genesis_command.split()
genesis_command_last = genesis_command_split[-1]
if not os.path.exists(os.path.expanduser(genesis_command_last)):
    print("Error: %s does not exist" % genesis_command_last)
    sys.exit(3)

if (is_fugaku):
    genesis_mpi_number = 8
else:
    is_mpi=False
    genesis_mpi_number = -1
    for options in genesis_command_split:
        if is_mpi:
            genesis_mpi_number = int(options)
            is_mpi=False
            break
        if options == "-n":
            is_mpi = True
        if options == "-np":
            is_mpi = True

    if int(genesis_mpi_number) < 0:
        genesis_mpi_number = genesis_command_split[-2]

    if not int(genesis_mpi_number) == 8:
        print("Error: Number of MPI processes %d should be 8" % int(genesis_mpi_number))
        sys.exit(3)

# if given path is relpath, change it to abspath
genesis_command_split[-1] = os.path.abspath(os.path.expanduser(genesis_command_last))
genesis_command = " ".join(genesis_command_split)

if genesis_command[-5:] == "atdyn":
    is_atdyn = True
elif genesis_command[-5:] == "spdyn":
    is_spdyn = True

###### setup test directories

if is_atdyn:
    test_dirs = getdirs(os.path.dirname(__file__) + "/test_gamd_atdyn")
elif is_spdyn:
    test_dirs = getdirs(os.path.dirname(__file__) + "/test_gamd_spdyn")

for dir in test_dirs:
    if not os.path.exists(dir):
        print("Error: %s, this test directory does not exist" % dir)
        sys.exit(3)

###### run tests
if (is_atdyn or is_spdyn):
    print("=======================================================================")
    print(" Regression tests for GaMD")
    print("=======================================================================")
    
    cwdname = os.getcwd()
    for test_each in test_dirs:
        itried = itried + 1
        os.chdir(cwdname)
        dirname = test_each
        if not os.path.isdir(dirname) : 
            continue
        os.chdir(dirname)
        
        # run MD
        print("-----------------------------------------------------------------------")
        print("Running %s..." % (dirname + "/"))

        for fl in glob.glob("test*"):
            os.remove(fl)
        if os.path.exists("out.gamd"):
            os.remove("out.gamd")
         
        inputname = "inp"
        testname = "test"
        if (is_fugaku):
            commandline = '%s sh -c \"%s %s 1> %s 2> error\"' % (mpiexec_command, genesis_path, inputname, testname)
        else:
            commandline = '%s %s 1> %s 2> error' % (genesis_command, inputname, testname)
        print("$ %s" % commandline)
        status = subprocess.getstatusoutput(commandline)
        
        if (status[0] > 0) and (status[0] != 1024):
            print()
            print("Aborted...")
            print()
            iaborted = iaborted + 1
            continue
        
        # parse the result
        #if status[0] == 0:
        print()
        print("Checking %s" % test_each)
        print()

        test = Genesis()
        is_empty = os.stat(testname).st_size == 0
        if (is_empty):
            print()
            print("Aborted...")
            print()
            iaborted = iaborted + 1
            continue

        test.read(testname)
        refname = "ref"
        tolerance_cur = tolerance
        if test.is_single and is_spdyn:
            tolerance_cur = tolerance_single
        if test.is_fujitsu:
            tolerance_cur = tolerance_fujitsu
        if test.is_fugaku and not test.is_single:
            tolerance_cur = tolerance_fujitsu

        tolerance_cur_virial = tolerance_cur*virial_ratio

        ref = Genesis()
        ref.read(refname)
 
        # check the result
        print()
        print("Checking %s" % test_each)
        print("Checking diff between %s and %s..." % (refname, testname))
        print()
        ref.test_diff(test, tolerance_cur, tolerance_cur_virial)
        print()
        if ref.is_passed:
            ipassed = ipassed + 1
        else:
            ifailed = ifailed + 1
        # post-cleaning
        #os.remove(testname)

###### finalization
if (itried > 0):
    print("-----------------------------------------------------------------------")
    print("Passed  %d / %d" % (ipassed,  itried))
    print("Failed  %d / %d" % (ifailed,  itried))
    print("Aborted %d / %d" % (iaborted, itried))
    print("-----------------------------------------------------------------------")

if iaborted > 0:
    sys.exit(2)
elif ifailed > 0:
    sys.exit(1)
else:
    sys.exit(0)

