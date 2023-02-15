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
import math
from genesis import *

############### DEFINITION ##################################
class Minfo(object):
    def __init__(self):  # initialization
        self.is_passed = False
        self.freq_data = []
        self.vec_data = []

    def read(self, filename):
        fid = open(filename, 'r')
        text_list = fid.readlines()
        fid.close()
        self.parse(text_list)

    def parse(self, text_list): 
        patternVIB  = re.compile('Vibrational Data')
        patternFREQ = re.compile('Vibrational Frequency')

        num_data_per_line = 5

        x = 0
        for x in range(0, len(text_list)):
            line = text_list[x]
            result = patternVIB.search(line)
            if result is not None:
                break

        x+=1
        line = text_list[x]
        #print str(x)+": "+line
        num_domain = 0
        result = re.search("Domain", line)
        if result:
            x+=1
            line = text_list[x]
            num_domain = int(text_list[x])

        nd = 0
        while nd < num_domain:
            xx = x
            for xx in range(x, len(text_list)):
               line = text_list[xx]
               result = patternFREQ.search(line)
               if result:
                 break

            x = xx + 1
            ndata = int(text_list[x])
            nline = ndata // num_data_per_line

            if ndata % num_data_per_line != 0:
               nline += 1 

            #print "ndata = %d" % ndata
            #print "nline = %d" % nline

            nl = 0
            while nl < nline:
               x += 1
               data_line = text_list[x].split(",")
               for i in data_line:
                  #print str(nl)+" "+i
                  self.freq_data.append((float(i)))
               nl += 1

            x += 1
            nm = 0
            while nm < ndata:
               #print "mode = %d" % nm
               x = x + 2
               nl = 0
               tmp = []
               while nl < nline:
                  x += 1
                  data_line = text_list[x].split(",")
                  for i in data_line:
                     #print str(nl)+" "+i
                     tmp = tmp + [(float(i))]
                  nl += 1
               self.vec_data.append(tmp)
               nm += 1

            nd += 1

    def test_diff(self, obj, tolerance): # compare energies
        is_failure = False
        mode_failure = []
        vec_failure  = []

        for imode in range(0, len(self.freq_data)):
            d = abs(self.freq_data[imode] - obj.freq_data[imode])
            if abs(self.freq_data[imode]) < 1e4: #min log is 1e-4
                ratio=d
            else:
                ebase=max(abs(self.freq_data[imode]),1.0)
                ratio = d/ebase
            tolerance2 = tolerance

            if abs(self.freq_data[imode]) < 1e4:
                tolerance2 = tolerance2*1e4
            if ratio > tolerance2:
                is_failure = True
                mode_failure.append(imode)

            for ivec in range(0, len(self.vec_data[imode])):
                d = abs(self.vec_data[imode][ivec]) - abs(obj.vec_data[imode][ivec])
                
                if abs(self.vec_data[imode][ivec]) < 1e4: #min log is 1e-4
                    ratio=d
                else:
                    ebase=max(abs(self.vec_data[imode][ivec]),1.0)
                    ratio = d/ebase

                if ratio > tolerance:
                    is_failure = True
                    vec_failure.append([imode, ivec])

        if is_failure:
            self.is_passed = False

            if len(mode_failure) > 0:
                print("Failure in frequency (tolerance = %4.2e)" % tolerance)
                for count in range(0, len(mode_failure)):
                    imode = mode_failure[count]
                    print("mode %d" % (imode+1))

                    sys.stdout.write("< ")
                    sys.stdout.write("%14s" % str(self.freq_data[imode]).rjust(14))
                    sys.stdout.write("\n")
                    
                    sys.stdout.write("> ")
                    sys.stdout.write("%14s" % str(obj.freq_data[imode]).rjust(14))
                    sys.stdout.write("\n\n")

            if len(vec_failure) > 0:
                print("Failure in vector (tolerance = %4.2e)" % tolerance)
                for count in range(0, len(vec_failure)):
                    imode = vec_failure[count][0]
                    ivec  = vec_failure[count][1]
                    print("mode %d" % (imode+1))
                    print("element %d" % (ivec+1))

                    sys.stdout.write("< ")
                    sys.stdout.write("%14s" % str(self.vec_data[imode][ivec]).rjust(14))
                    sys.stdout.write("\n")
                    
                    sys.stdout.write("> ")
                    sys.stdout.write("%14s" % str(obj.vec_data[imode][ivec]).rjust(14))
                    sys.stdout.write("\n\n")

        else:
            self.is_passed = True
            print("Passed (tolerance = %4.2e(freq, vec))" % (tolerance))

############### DEFINITION ##################################
def getdirs(path):
    test_dirs = []
    for test_dir in os.listdir(path):
        test_dir_path = os.path.join(path,test_dir)
        #print "%s" % test_dir_path
        if os.path.isdir(test_dir_path):
            if os.path.exists(test_dir_path + "/inp") and os.path.exists(test_dir_path + "/ref1"):
                test_dirs.append(test_dir_path)
    test_dirs.sort()
    return test_dirs

############### MAIN ########################################

###### initialization

#os.environ["OMP_NUM_THREADS"] = "1"
tolerance = 1.0e-8 # relative energy difference (diff/abs(e))
tolerance_single = 3.0e-5 # relative energy difference (diff/abs(e))
virial_ratio = 1.0e2

tolerance_virial = tolerance*virial_ratio 

ipassed = 0
ifailed = 0
iaborted = 0
itried = 0
num = 0

is_atdyn = False
is_spdyn = False
is_parallelio = False
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
    $ ./test_vib.py ["genesis command"] [parallel_io or tolerance_value or directories]
    
    examples:
    # run tests using the default command ("mpirun -np 8 atdyn")
    $ ./test_vib.py

    # run atdyn tests
    $ ./test_vib.py "mpirun -np 8 /path/to/genesis/bin/atdyn"

    # run with specfic tolerance value
    $ ./test_vib.py "mpirun -np 8 /path/to/genesis/bin/atdyn" 0.1

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
    genesis_mpi_number = genesis_command.split()[-2]
    if not int(genesis_mpi_number)%4 == 0:
        print("Error: %d should be multiplier of 4" % int(genesis_mpi_number))
        sys.exit(3)

# if given path is relpath, change it to abspath
genesis_command_split[-1] = os.path.abspath(os.path.expanduser(genesis_command_last))
genesis_command = " ".join(genesis_command_split)

if genesis_command[-5:] == "atdyn":
    is_atdyn = True
elif genesis_command[-5:] == "spdyn":
    is_spdyn = True

if is_spdyn:
    print("Error: test_vib for spdyn is not available")
    sys.exit(3)
  
###### setup test directories

if is_atdyn:
    test_dirs = getdirs(os.path.dirname(os.path.abspath(__file__)) + "/test_vib")

for dir in test_dirs:
    if not os.path.exists(dir):
        print("Error: %s, this test directory does not exist" % dir)
        sys.exit(3)

###### run tests
if (is_atdyn or is_spdyn):
    print("=======================================================================")
    print(" Regression tests for VIBRATION")
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

        minfo="vib.minfo"
        if os.path.exists(minfo):
          os.remove(minfo)

        minfodir="minfo.files"
        if os.path.exists(minfodir):
          for fl in glob.glob(minfodir + "/*.minfo"):
              os.remove(fl)
        
        inputname = "inp"
        outname = "log"
        if (is_fugaku):
            commandline = '%s sh -c \"%s %s 1> %s 2> error\"' % (mpiexec_command, genesis_path, inputname, outname)
        else:
            commandline = '%s %s 1> %s 2> error' % (genesis_command, inputname, outname)
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

        num = 1
        ipassed_replica = 0
        ifailed_replica = 0
        while num < 9:
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
            ref.test_diff(test, tolerance, tolerance_virial)
            print()

            if ref.is_passed:
                ipassed_replica = ipassed_replica + 1
            else:
                ifailed_replica = ifailed_replica + 1

            num += 1
            # post-cleaning
            #os.remove(testname)
	    
        ref = Minfo()
        refname = "ref.minfo"
        ref.read(refname)

        test = Minfo()
        testname = "vib.minfo"
        test.read(testname)

        # check the result
        print()
        print("Checking diff between %s and %s..." % (refname, testname))
        print()

        ref.test_diff(test, tolerance)

        if (ifailed_replica == 0 and ref.is_passed):
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

