#!/usr/bin/python
# coding: utf-8

#
# A parser script for GENESIS output style
#
# (c) Copyright 2022 RIKEN. All rights reserved.
#

import re
import sys
import copy
import math

############### DEFINITION ####################################################
class Genesis(object):
    def __init__(self):  # initialization
        # public attributes
        self.dict_text = {}
        self.dict_data = {}
        self.dict_error = {}
        self.is_passed = False
        self.is_fujitsu = False
        self.is_fugaku = False
        self.is_single = False
        self.is_gpu    = False

    def delete_last(self):
        for key in list(self.dict_data.keys()):
            self.dict_text[key].pop()
            self.dict_data[key].pop()

    def delete_first(self):
        for key in list(self.dict_data.keys()):
            self.dict_text[key].pop(0)
            self.dict_data[key].pop(0)

    def read(self, filename):
        fid = open(filename, 'r')
        text_list = fid.readlines()
        fid.close()
        self.parse(text_list)

    def parse(self, text_list): 
        # parse titles
        title = []
        text = []
        data = []
        data_each = []
        #patternMD = re.compile('^INFO:')
        patternMD = re.compile('^INFO.')
        patternMINUS = re.compile('-')
        patternNONBOND = re.compile('^  nonbonding')
        patternCPU = re.compile('^  cpu model')
        patternFUJITSU = re.compile('SPARC')
        patternFUGAKU = re.compile('Fugaku')
        patternPRECISION = re.compile('^  precision')
        patternSINGLE = re.compile('single')
        patternMIXED = re.compile('mixed')
        patternGPU = re.compile('GPU')
        patternSTEP0 = re.compile('^\[STEP0')
        patternSTEP1 = re.compile('^\[STEP1')
        is_header = False
        is_md = False
        for line in text_list:
            if is_md:
                result = patternMD.search(line)
                if result is None:
                    continue
                else:  # patternMD found
                    line_sub = patternMD.sub('', line.rstrip('\n'))
                    line_sub = patternMINUS.sub(' -', line_sub)
                    data_each = line_sub.split()
                    text.append(data_each)
                    tmp = []
                    for i in range(len(data_each)):
                        result = False
                        if len(data_each[i]) >= 2 :
                            if data_each[i][-1].isdigit() & data_each[i][0].isalpha() :
                                result = True

                        if result is True : # number ?
                            tmp = tmp + [(float(999999))] 
                        elif data_each[i].replace('.','',1).isdigit(): # positive number
#                            print "test-posi %s" % data_each[i]
                            tmp = tmp + [(float(data_each[i]))]
                        elif data_each[i][0]=="-" and data_each[i][1:].replace('.','',1).isdigit(): # negative number 
#                            print "test-nega %s" % data_each[i]
                            tmp = tmp + [(float(data_each[i]))]
                        elif data_each[i][-1].isdigit() and data_each[i][0].isdigit(): # number ?
#                            print "test-nani %s" % data_each[i]
                            tmp = tmp + [(float(data_each[i]))]
                        else: # Not a Number !
                            tmp = tmp + [(float(999999))] 
                    data.append(tmp)
            else:
                result = patternMD.search(line)
                if result is not None:
                    line_sub = patternMD.sub('', line.rstrip('\n'))
                    line_split = line_sub.split()
                    title = title + line_split
                    is_md = True

            if is_header:
                result = patternSTEP1.search(line)
                if result is not None:
                    is_header = False
                    continue
                result = patternCPU.search(line)
                if result is not None:
                    self.is_fujitsu = patternFUJITSU.search(line)
                if result is not None:
                    self.is_fugaku = patternFUGAKU.search(line)
                result = patternPRECISION.search(line)
                if result is not None:
                    self.is_single = patternSINGLE.search(line)
                    if not self.is_single:
                        self.is_single = patternMIXED.search(line)
                result = patternNONBOND.search(line)
                if result is not None:
                    self.is_gpu = patternGPU.search(line)
            else:
                result = patternSTEP0.search(line)
                if result is not None:
                    is_header = True

        # append to the dictionary
        self.dict_append(title, text, data)

    def dict_append(self, title, text, data):
        if len(self.dict_text) == 0:
            self.dict_text = dict.fromkeys(title, [])
        if len(self.dict_data) == 0:
            self.dict_data = dict.fromkeys(title, [])
        text_transpose = list(map(list, list(zip(*text))))
        data_transpose = list(map(list, list(zip(*data))))
        for i in range(len(title)):
            self.dict_text[title[i]] = self.dict_text[title[i]] + text_transpose[i]
            self.dict_data[title[i]] = self.dict_data[title[i]] + data_transpose[i]

    def test_diff(self, obj, tolerance, tolerance_virial): # compare energies
        # test MD steps
        keys = set(self.dict_data.keys()) & set(obj.dict_data.keys())
        is_failure = False
        dict_failure = dict.fromkeys(keys, False)
        nstep_failure = 0
        patternPRESS = re.compile('PRESS')

        nstep = len(self.dict_data['STEP'])
        for istep in range(nstep):
            for key in keys:
                d = abs(self.dict_data[key][istep] - obj.dict_data[key][istep])
                if abs(self.dict_data[key][istep]) < 1e4: #min log is 1e-4
                    ratio=d
                else:
                    ebase=max(abs(self.dict_data[key][istep]),1.0)
                    ratio = d/ebase
                if (key == "VIRIAL"):
                    tolerance2 = tolerance_virial
                elif (patternPRESS.search(key)):
                    tolerance2 = tolerance_virial
                else:
                    tolerance2 = tolerance
                if abs(self.dict_data[key][istep]) < 1e4:
                    tolerance2 = tolerance2*1e4
                if ratio > tolerance2:
                    is_failure = True
                    dict_failure[key] = True
                    nstep_failure = istep
            if is_failure:
                break

        if is_failure:
            self.is_passed = False
            print("Failure at step %d (tolerance = %4.2e(ene), %4.2e(virial))" % (self.dict_data['STEP'][nstep_failure], tolerance,tolerance_virial))
            nstep_max = min([nstep_failure + 3, nstep])
            for istep in range(nstep_failure, nstep_max):
                print("Step %d" % (self.dict_data['STEP'][istep]))

                sys.stdout.write("  ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % key.rjust(14))
                sys.stdout.write("\n")

                sys.stdout.write("< ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % self.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n")
                
                sys.stdout.write("> ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % obj.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n\n")
        else:
            self.is_passed = True
            print("Passed (tolerance = %4.2e(ene), %4.2e(virial))" % (tolerance, tolerance_virial))

    def test_diff_TMD(self, obj, tolerance, tolerance_virial): # compare energies for TMD
        # test MD steps
        keys = set(self.dict_data.keys()) & set(obj.dict_data.keys())
        is_failure = False
        dict_failure = dict.fromkeys(keys, False)
        nstep_failure = 0
        patternPRESS = re.compile('PRESS')

        nstep = len(self.dict_data['STEP'])
        for istep in range(nstep):
            for key in keys:
                d = abs(self.dict_data[key][istep] - obj.dict_data[key][istep])
                if abs(self.dict_data[key][istep]) < 1e4: #min log is 1e-4
                    ratio=d
                else:
                    ebase=max(abs(self.dict_data[key][istep]),1.0)
                    ratio = d/ebase
                if (key == "RESTRAINT_TOTAL"):
                    continue
                if (key == "TOTAL_ENE"):
                    continue
                if (key == "POTENTIAL_ENE"):
                    continue
                elif (key == "VIRIAL"):
                    tolerance2 = tolerance_virial
                elif (patternPRESS.search(key)):
                    tolerance2 = tolerance_virial
                else:
                    tolerance2 = tolerance
                if abs(self.dict_data[key][istep]) < 1e4:
                    tolerance2 = tolerance2*1e4
                if ratio > tolerance2:
                    is_failure = True
                    dict_failure[key] = True
                    nstep_failure = istep
            if is_failure:
                break

        if is_failure:
            self.is_passed = False
            print("Failure at step %d (tolerance = %4.2e(ene), %4.2e(virial))" % (self.dict_data['STEP'][nstep_failure], tolerance,tolerance_virial))
            nstep_max = min([nstep_failure + 3, nstep])
            for istep in range(nstep_failure, nstep_max):
                print("Step %d" % (self.dict_data['STEP'][istep]))

                sys.stdout.write("  ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % key.rjust(14))
                sys.stdout.write("\n")

                sys.stdout.write("< ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % self.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n")
                
                sys.stdout.write("> ")
                for key in keys:
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % obj.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n\n")
        else:
            self.is_passed = True
            print("Passed (tolerance = %4.2e(ene), %4.2e(virial))" % (tolerance, tolerance_virial))

############### TEST ##########################################################
if __name__ == '__main__':
    text = '''Update_Pairlist_Pbc> Memory for Pairlist was allocated
Update_Pairlist_Pbc> Pairlist was Updated
[STEP4] Compute Single Point Energy for Molecules
 
            STEP            BOND           ANGLE        DIHEDRAL        IMPROPER
         VDWAALS           ELECT    UREY-BRADLEY            CMAP       RESTRAINT
 --------------- --------------- --------------- --------------- ---------------
               0       3590.8442       1851.3287          0.0000          0.0000
       4565.8064     -39052.9065          0.0000          0.0000          0.0000

[STEP5] Perform Molecular Dynamics Simulation
 
Update_Pairlist_Pbc> Pairlist was Updated
INFO:       STEP            TIME       TOTAL_ENE   POTENTIAL_ENE     KINETIC_ENE            RMSG            BOND           ANGLE        DIHEDRAL        IMPROPER         VDWAALS           ELECT    UREY-BRADLEY            CMAP       RESTRAINT            BOXX            BOXY            BOXZ          VOLUME     TEMPERATURE         PRESSXX         PRESSYY         PRESSZZ          VIRIAL        PRESSURE
 --------------- --------------- --------------- --------------- ---------------
INFO:          1          0.0010     -20366.0251     -29044.9272       8678.9021         19.3525       3590.8442       1851.3287          0.0000          0.0000       4565.8064     -39052.9065          0.0000          0.0000          0.0000         45.8363         45.8363         45.8363      96300.5261        315.9647        -58.0979       -470.2127        939.1533      -5593.6679        136.9476

Update_Pairlist_Pbc> Memory for Pairlist was allocated
Update_Pairlist_Pbc> Pairlist was Updated
INFO:          2          0.0020     -20357.3017     -28943.4688       8586.1671         19.7050       3666.5101       1848.7082          0.0000          0.0000       4570.7423     -39029.4294          0.0000          0.0000          0.0000         45.8363         45.8363         45.8363      96300.5261        312.5886        187.0513       -715.0704        842.0369      -5577.1568        104.6726

[STEP6] Deallocate Arrays
 
Output_Time> Timer profile of each rank
   Rank       Ebond      Enbond       Integ        List       Total
      0       0.001       0.172       0.001       1.021       1.706
      1       0.000       0.176       0.001       1.017       1.708
      2       0.000       0.178       0.001       1.015       1.709
      3       0.000       0.173       0.001       1.020       1.709
      4       0.000       0.178       0.001       1.015       1.708
      5       0.002       0.175       0.001       1.018       1.710
      6       0.001       0.180       0.001       1.012       1.708
      7       0.000       0.174       0.001       1.019       1.708

Output_Time> Averaged timer profile
  total time      =       1.708
    setup         =       0.976
    dynamics      =       0.732
      energy      =       0.179
      integrator  =       0.001
      pairlist    =       1.017
  energy           
    bond          =       0.000
    angle         =       0.000
    dihedral      =       0.000
    nonbond       =       0.176
      pme real    =       0.146
      pme recip   =       0.029
    restraint     =       0.000
  integrator       
    constraint    =       0.000
    update        =       0.000
    comm1         =       0.000
    comm2         =       0.000
'''
    ene = Genesis()

    # test parse
    ene.parse(text.split("\n"))
    assert len(ene.dict_data) == 25
    assert ene.dict_text['ELECT'][0] == "-39052.9065"
    assert ((ene.dict_data['ELECT'][0] + 39052.9065) < 10**(-3))

    # test diff
    ene2 = copy.deepcopy(ene)
    ene2.dict_text['ELECT'][0] = "-39052.0293"
    ene2.dict_data['ELECT'][0] = -39052.0293
    ene.test_diff(ene2, 0.0001, 0.001)
    ene.delete_last()


