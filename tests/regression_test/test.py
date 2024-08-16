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
    for system_dir in os.listdir(path):
        system_dir_path = os.path.join(path,system_dir)
        if os.path.isdir(system_dir_path):
            for test_dir in os.listdir(system_dir_path):
                test_dir_path = os.path.join(system_dir_path,test_dir)
                if os.path.isdir(test_dir_path):
                    if os.path.exists(test_dir_path + "/inp") and os.path.exists(test_dir_path + "/ref"):
                        test_dirs.append(test_dir_path)
                    if os.path.exists(test_dir_path + "/pio.inp"):
                        test_dirs.append(test_dir_path)
    test_dirs.sort()
    return test_dirs

############### MAIN ########################################

###### initialization

tolerance = 1.0e-6 # relative energy difference (diff/abs(e))
tolerance_atdyn = 1.0e-8 # relative energy difference (diff/abs(e))
tolerance_tmd = 5.0e-6 # relative energy difference (diff/abs(e))
tolerance_fujitsu = 5.0e-6 # relative energy difference (diff/abs(e))
tolerance_gpu_respa = 5.0e-6 # relative energy difference (diff/abs(e))
tolerance_fujitsu_weak = 1.0e-5 # relative energy difference (diff/abs(e))
tolerance_press_rmsd   = 3.0e-5 # relative energy difference (diff/abs(e))
tolerance_single = 3.0e-5 # relative energy difference (diff/abs(e))
virial_ratio = 2.0e2
#os.environ["OMP_NUM_THREADS"] = "1"

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
    genesis_command = 'mpirun -np 8 atdyn'
elif len(sys.argv) == 2:
    genesis_command = sys.argv[1]
else:
    genesis_command = sys.argv[1]
    is_number = re.compile(r'^[+-]?(\d*\.\d+|\d+\.?\d*)([eE][+-]?\d+|)\Z')
    if is_number.match(sys.argv[2]):
        tolerance = float(sys.argv[2])
    elif sys.argv[2] == "parallel_io":
        is_parallelio = True
    elif sys.argv[2] == "gpu":
        is_gpu = True
        if len(sys.argv) == 4:
            if is_number.match(sys.argv[3]):
                tolerance = float(sys.argv[3])
    elif sys.argv[2] == "fugaku":
        mpiexec_command = genesis_command.split(' ',2)[0]
        genesis_path = genesis_command.split(' ',2)[1]
        is_fugaku = True
        if len(sys.argv) == 4:
            if is_number.match(sys.argv[3]):
                tolerance = float(sys.argv[3])
    else:
        test_dirs = sys.argv[2:]

tolerance_virial = tolerance*virial_ratio 

if (genesis_command == "-h") or (genesis_command == "--help"):
    print("""    usage:
    $ ./test.py ["genesis command"] [parallel_io or tolerance_value or directories]
    
    examples:
    # run tests using the default command ("mpirun -np 8 atdyn")
    $ ./test.py

    # run atdyn tests
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/atdyn"

    # run spdyn tests
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/spdyn"

    # run parallel IO tests
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/spdyn" parallel_io

    # run gpu tests
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/spdyn" gpu

    # run with specfic tolerance value
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/atdyn" 0.1

    # run specfic directories
    $ ./test.py "mpirun -np 8 /path/to/genesis/bin/atdyn" test_common/jac_param27/LEAP_CUTOFF
    """)
    sys.exit(3)

genesis_command_split = genesis_command.split()
genesis_command_last = genesis_command_split[-1]
if not os.path.exists(os.path.expanduser(genesis_command_last)):
    print("Error: %s does not exist" % genesis_command_last)
    sys.exit(3)

# if given path is relpath, change it to abspath
genesis_command_split[-1] = os.path.abspath(os.path.expanduser(genesis_command_last))
genesis_command = " ".join(genesis_command_split)

if genesis_command[-5:] == "atdyn":
    is_atdyn = True
elif genesis_command[-5:] == "spdyn":
    is_spdyn = True

###### setup test directories

if len(test_dirs) == 0:
    if is_parallelio and is_spdyn:
        test_dirs = getdirs(os.path.dirname(os.path.abspath(__file__)) + "/test_parallel_IO")
    elif is_parallelio and (not is_spdyn):
        print("Error: parallel_io is supported by only spdyn, please specify spdyn in command line")
        sys.exit(3)
    elif is_spdyn:
        test_dirs =  getdirs(os.path.dirname(os.path.abspath(__file__)) + "/test_spdyn")
    elif is_atdyn:
        test_dirs =  getdirs(os.path.dirname(os.path.abspath(__file__)) + "/test_atdyn")

for dir in test_dirs:
    if not os.path.exists(dir):
        print("Error: %s, this test directory does not exist" % dir)
        sys.exit(3)

###### check the number of MPI processors
if (len(sys.argv) == 2 and is_atdyn):
    genesis_mpi_number = 1

elif (is_fugaku):
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

    if (is_parallelio and is_spdyn):
        if int(genesis_mpi_number) != 8:
            print("Error: Number of MPI processes should be 8 for parallel_io tests")
            sys.exit(3)
    else:
        if not int(genesis_mpi_number) in [1, 2, 4, 8]:
            print("Error: Number of MPI processes should be 1, 2, 4, or 8")
            sys.exit(3)

###### run tests
if (is_atdyn or is_spdyn) and (not is_parallelio):
    print("=======================================================================")
    print(" Regression tests for MD")
    print("=======================================================================")
    
    cwdname = os.getcwd()
    for test_each in test_dirs:
        os.chdir(cwdname)
        dirname = test_each
        resultvver = re.compile(r'VVER_[A-Za-z_]+([1248])').search(dirname);
        resultvres = re.compile(r'VRES_[A-Za-z_]+([1248])').search(dirname);
        if ("WATER" in dirname) and (is_gpu) : 
            continue
        if ("CUTOFF" in dirname) and (is_gpu) : 
            continue
        if ("Go" in dirname) and (is_gpu) : 
            continue
        if ("All_atom" in dirname) and (is_gpu) : 
            continue
        if ("cago" in dirname) and (is_gpu) : 
            continue
        if (resultvver and (is_spdyn)):
            if (int(genesis_mpi_number) != int(resultvver.group(1))) :
                continue
        if (resultvres and (is_spdyn)):
            if (int(genesis_mpi_number) != int(resultvres.group(1))) :
                continue

        itried = itried + 1
        if not os.path.isdir(dirname) : 
            continue
        os.chdir(dirname)
        
        # run MD
        print("-----------------------------------------------------------------------")
        print("Running %s..." % (dirname + "/"))
        
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
        if is_atdyn:
            tolerance_cur = tolerance_atdyn
        if is_atdyn and ("TMD" in dirname):
            tolerance_cur = tolerance_tmd
        if is_spdyn and ("TMD" in dirname) and not test.is_single:
            tolerance_cur = tolerance_tmd
        if not test.is_single and is_spdyn and is_gpu and ("VRES" in dirname): 
            tolerance_cur = tolerance_gpu_respa
        if test.is_fujitsu:
            tolerance_cur = tolerance_fujitsu
        if is_fugaku and not test.is_single:
            tolerance_cur = tolerance_fujitsu
        if test.is_fujitsu and is_spdyn and ("tip4" in dirname):
            tolerance_cur = tolerance_fujitsu_weak
        if is_fugaku and is_spdyn and ("tip4" in dirname):
            tolerance_cur = tolerance_fujitsu_weak
        if is_atdyn and ("REST-RMSD_EXT" in dirname):
            tolerance_cur = tolerance_press_rmsd

        tolerance_cur_virial = tolerance_cur*virial_ratio

        ref = Genesis()
        ref.read(refname)
        
        # check the result
        print()
        print("Checking %s" % test_each)
        print("Checking diff between %s and %s..." % (refname, testname))
        print()
        if ("TMD" in dirname and test.is_single):
            ref.test_diff_TMD(test, tolerance_cur, tolerance_cur_virial)
        else:
            ref.test_diff(test, tolerance_cur, tolerance_cur_virial)

        print()
        if ref.is_passed:
            ipassed = ipassed + 1
        else:
            ifailed = ifailed + 1
        # post-cleaning
        #os.remove(testname)

###### run parallel i/o tests
if is_parallelio and is_spdyn:
    print("=======================================================================")
    print("Parallel I/O Tests")
    print("=======================================================================")


    genesis_command_prst_setup_split = genesis_command_split
    genesis_prst_setup=genesis_command_prst_setup_split[-1].replace('spdyn','prst_setup');
    genesis_prst_setup = os.path.abspath(os.path.expanduser(genesis_prst_setup));
    if not os.path.exists(os.path.expanduser(genesis_prst_setup)):
        print("Error: %s does not exist" % genesis_prst_setup)
        sys.exit(3)

    genesis_command_prst_setup_split[-1] = genesis_prst_setup
    genesis_prst_setup = " ".join(genesis_command_prst_setup_split)
        
    cwdname = os.getcwd()
    
    for test_each in test_dirs:
        dirname = test_each
        os.chdir(dirname)
        itried = itried + 1
        for fl in glob.glob("dd*.rst"):
            os.remove(fl)
        if os.path.isdir("./cache"): 
            shutil.rmtree("cache")
        os.mkdir("./cache")

        # run MD
        print("-----------------------------------------------------------------------")
        print("Building Parallel I/O rst %s..." % (dirname + "/"))
        commandline = '%s prst.inp 1> log 2> error_prst' % genesis_prst_setup
        status1 = subprocess.getstatusoutput(commandline)
        if (status1[0] > 0) and (status1[0] != 1024):
            print()
            print("Aborted...")
            print()
            iaborted = iaborted + 1
            os.chdir(cwdname)
            continue
        
        print("Done ...")
        print("Running single version %s..." % (dirname + "/"))
        commandline = '%s single.inp 1> ref 2> error_single' % genesis_command
        status2 = subprocess.getstatusoutput(commandline)
        if (status2[0] > 0) and (status2[0] != 1024):
            print()
            print("Aborted...")
            print()
            iaborted = iaborted + 1
            os.chdir(cwdname)
            continue

        print("Done ...")
        testname = "test"
        commandline = '%s pio.inp 1> %s 2> error' % (genesis_command, testname)
        print("$ %s" % commandline)
        status4 = subprocess.getstatusoutput(commandline)
        
        if (status4[0] > 0) and (status4[0] != 1024):
            print()
            print("Aborted...")
            print()
            iaborted = iaborted + 1
            os.chdir(cwdname)
            continue
        
        # parse the result
        #if status[0] == 0:
        refname = "ref"
        ref = Genesis()
        ref.read(refname)
        test = Genesis()
        test.read(testname)
        
        # check the result
        print()
        print("Checking diff between %s and %s..." % (refname, testname))
        print()
        ref.test_diff(test, tolerance, tolerance_virial)
        print()
        if ref.is_passed:
            ipassed = ipassed + 1
        else:
            ifailed = ifailed + 1
        # post-cleaning
        #os.remove(testname)
        os.chdir(cwdname)

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

