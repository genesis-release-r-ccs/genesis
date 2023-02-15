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
from fep import *


############### DEFINITION ##################################
def getdirs(path):
    test_dirs = []
    for forcefield_dir in os.listdir(path):
        forcefield_dir_path = os.path.join(path,forcefield_dir)
        if os.path.isdir(forcefield_dir_path):
            for test_dir in os.listdir(forcefield_dir_path):
                test_dir_path = os.path.join(forcefield_dir_path,test_dir)
                if os.path.isdir(test_dir_path):
                    if os.path.exists(test_dir_path + "/inp") and (os.path.exists(test_dir_path + "/ref") or os.path.exists(test_dir_path + "/ref1")):
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
tolerance_single = 1.0e-4 # relative energy difference (diff/abs(e))
virial_ratio = 1.0e2

ipassed = 0
ifailed = 0
iaborted = 0
itried = 0

is_atdyn = False
is_spdyn = False
is_parallelio = False
is_gpu = False
is_fugaku = False

test_dirs = []

###### parse command line

if len(sys.argv) == 1:
    genesis_command = 'mpirun -np 8 spdyn'
elif len(sys.argv) == 2:
    genesis_command = sys.argv[1]
else:
    genesis_command = sys.argv[1]
    is_number = re.compile('^[+-]?(\d*\.\d+|\d+\.?\d*)([eE][+-]?\d+|)\Z')
    if is_number.match(sys.argv[2]):
        tolerance = float(sys.argv[2])
    elif sys.argv[2] == "gpu":
        is_gpu = True
#         if len(sys.argv) == 4:
#             if is_number.match(sys.argv[3]):
#                 tolerance = float(sys.argv[3])
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
    # run tests using the default command ("mpirun -np 8 spdyn")
    $ ./test_fep.py

    # run spdyn tests
    $ ./test_fep.py "mpirun -np 8 /path/to/genesis/bin/spdyn"

    # run gpu tests
    $ ./test_fep.py "mpirun -np 8 /path/to/genesis/bin/spdyn" gpu

    # run with specfic tolerance value
    $ ./test_fep.py "mpirun -np 8 /path/to/genesis/bin/spdyn" 0.1

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

if genesis_command[-5:] == "spdyn":
    is_spdyn = True
if genesis_command[-5:] == "atdyn":
    print("FEP is not available in atdyn. Use spdyn")
    sys.exit(3)

###### setup test directories

if is_spdyn:
    test_dirs = getdirs(os.path.dirname(os.path.abspath(__file__)) + "/test_fep")

for dir in test_dirs:
    if not os.path.exists(dir):
        print("Error: %s, this test directory does not exist" % dir)
        sys.exit(3)

###### run tests
if (is_spdyn):
    print("=======================================================================")
    print(" Regression tests for FEP")
    print("=======================================================================")
    
    cwdname = os.getcwd()
    for test_each in test_dirs:
        dirname = test_each
        # skip CUTOFF and minimization if gpu is on
        if "CUTOFF" in dirname or "MIN" in dirname: 
            if is_gpu:
                continue
        if "REMD" in dirname:
            if is_fugaku:
                continue
        os.chdir(cwdname)
        if not os.path.isdir(dirname) : 
            continue
        os.chdir(dirname)
        itried = itried + 1

        # run MD
        print("-----------------------------------------------------------------------")
        print("Running %s..." % (dirname + "/"))

        for fl in glob.glob("test*"):
            os.remove(fl)
        for fl in glob.glob("out*.fepout"):
            os.remove(fl)
         
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

        if "REMD" in dirname or "FEPREST" in dirname: 
            num = 1
            ipassed_replica = 0
            ifailed_replica = 0

            while num < 5:
                ref = Genesis()
                refname = "ref%d" % num
                ref.read(refname)
                testname = "test%d" % num
                test = Genesis()
                test.read(testname)
            
                # check the result
                print()
                print("Checking diff between %s and %s..." % (refname, testname))
                print()
                #ref.test_diff_energies(test, tolerance)
                ref.test_diff(test, tolerance_cur, tolerance_cur_virial)
                print()

                if ref.is_passed:
                    ipassed_replica = ipassed_replica + 1
                else:
                    ifailed_replica = ifailed_replica + 1

                num += 1
                # post-cleaning
                #os.remove(testname)
                
            if (ifailed_replica == 0):
                ipassed = ipassed + 1
            else:
                ifailed = ifailed + 1

        else:
            ref = Genesis()
            ref.read(refname)
 
            # check the result
            print()
            print("Checking diff between %s and %s..." % (refname, testname))
            print()
            ref.test_diff(test, tolerance_cur, tolerance_cur_virial)
            print()

            if os.path.exists("ref.fepout"):
                test_fepout = Fep()
                testname = "out.fepout"
                test_fepout.read(testname)
                tolerance_cur = tolerance
                if test.is_single and is_spdyn:
                    tolerance_cur = tolerance_single

                refname = "ref.fepout"
                ref_fepout = Fep()
                ref_fepout.read(refname)

                # check the fepout result
                print("Checking diff between %s and %s..." % (refname, testname))
                print()
                ref_fepout.test_diff(test_fepout, tolerance_cur)
                print()
                if ref.is_passed and ref_fepout.is_passed:
                    ipassed = ipassed + 1
                else:
                    ifailed = ifailed + 1
            else:
                if ref.is_passed:
                    ipassed = ipassed + 1
                else:
                    ifailed = ifailed + 1


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

