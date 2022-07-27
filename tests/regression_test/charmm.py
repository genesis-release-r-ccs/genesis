#!/usr/bin/python
# coding: utf-8

#
# A parser script for CHARMM output style
#
# (c) Copyright 2022 RIKEN. All rights reserved.
#

import re
import sys
import copy

############### DEFINITION ####################################################
class Charmm(object):
    def __init__(self):  # initialization
        # public attributes
        self.dict_text = {}
        self.dict_data = {}
        self.is_passed = False

    def delete_last(self):
        for key in list(self.dict_data.keys()):
            self.dict_text[key].pop()
            self.dict_data[key].pop()

    def delete_first(self):
        for key in list(self.dict_data.keys()):
            self.dict_text[key].pop(0)
            self.dict_data[key].pop(0)

    def extract_first(self):
        for key in list(self.dict_data.keys()):
            self.dict_text[key] = [self.dict_text[key][0]]
            self.dict_data[key] = [self.dict_data[key][0]]

    def read(self, filename):
        fid = open(filename, 'r')
        text_list = fid.readlines()
        fid.close()
        self.parse(text_list)

    def parse(self, text_list): 
        # parse titles
        title = []
        patternMD  = re.compile(' CHARMM>') #for CHARMM
        patternMD2 = re.compile('Perform Molecular Dynamics Simulation|Perform Energy Minimization') #for GENESIS
        patternLABEL = re.compile('DYNA.*:')
        patternDYNA = re.compile('DYNA.*>')
        patternCONSTR = re.compile('DYNA CONSTR:')
        is_md = False
        for line in text_list:
            if is_md:
                result = patternDYNA.search(line)
                if result is not None:
                    break
                result = patternCONSTR.search(line)
                if result is not None:
                    continue
                result = patternLABEL.search(line)
                if result is None:
                    continue
                else: # patternLABEL found
                    line_sub = patternLABEL.sub('', line.rstrip('\n'))
                    line_split = line_sub.split()
                    title = title + line_split
            else:
                result = patternMD.search(line)
                if result is not None:
                    is_md = True
                result = patternMD2.search(line)
                if result is not None:
                    is_md = True
        
        # parse data
        text = []
        data = []
        data_each = []
        is_md = False
        patternMD  = re.compile(' CHARMM>') #for CHARMM
        patternMD2 = re.compile('Perform Molecular Dynamics Simulation|Perform Energy Minimization') #for GENESIS
        patternDYNA = re.compile('DYNA.*>')
        patternCONSTR = re.compile('DYNA CONSTR>')
        #patternDATA = re.compile('(\S+?\.\S{5}|\d+)')
        #patternDATA = re.compile('([\s-]+[\d\.]+)')
        patternDATA = re.compile('(-?\d+?\.\d{4,5}|\d+|NaN)')
        for line in text_list:
            if is_md:
                result = patternDYNA.search(line)
                if result is None:
                    continue
                else:  # patternDYNA found
                    result = patternCONSTR.search(line)
                    if result is not None:
                        continue
                    line_sub = patternDYNA.sub('', line.rstrip('\n'))
                    line_split = patternDATA.findall(line_sub)
                    data_each = data_each + line_split
                    # print line
                    # print data_each
                    # print
                    if len(data_each) == len(title):
                        text.append(data_each)
                        # data.append([float(item) for item in data_each])
                        tmp = []
                        for i in range(len(data_each)):
                            if data_each[i][-1].isdigit(): # number ?
                                tmp = tmp + [(float(data_each[i]))]
                            else: # Not a Number !
                                tmp = tmp + [(float(999999))] 
                        data.append(tmp)
                        data_each = []
            else:
                result = patternMD.search(line)
                if result is not None:
                    is_md = True
                result = patternMD2.search(line)
                if result is not None:
                    is_md = True

        # remove duplicate steps
        text2 = [];
        data2 = [];
        steps = [];
        for i in range(len(text)):
            if not data[i][0] in steps:
                steps.append(data[i][0])
                text2.append(text[i])
                data2.append(data[i])
            
        # append to the dictionary
        self.dict_append(title, text2, data2)

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

        # count IMAGE terms because these are not written in GENESIS output
        if 'IMNBvdw' in self.dict_data and 'VDWaals' in self.dict_data:
            nstep = len(self.dict_data['VDWaals'])
            for istep in range(nstep):
                self.dict_data['VDWaals'][istep] = self.dict_data['VDWaals'][istep] + self.dict_data['IMNBvdw'][istep]
                self.dict_text['VDWaals'][istep] = "%f" % self.dict_data['VDWaals'][istep]
            # del self.dict_data['IMNBvdw']
            # del self.dict_text['IMNBvdw']
        if 'IMELec' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['IMELec'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
        #     del self.dict_data['IMELec']
        #     del self.dict_text['IMELec']

        # if self.dict_data.has_key('IMHBnd'):
        #     del self.dict_data['IMHBnd']
        #     del self.dict_text['IMHBnd']
        # if self.dict_data.has_key('RXNField'):
        #     del self.dict_data['RXNField']
        #     del self.dict_text['RXNField']
        # if self.dict_data.has_key('EXTElec'):
        #     del self.dict_data['EXTElec']
        #     del self.dict_text['EXTElec']

        # count EWALD terms
        if 'EWKSum' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['EWKSum'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
        if 'EWSElf' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['EWSElf'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
        if 'EWEXcl' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['EWEXcl'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
        if 'EWQCor' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['EWQCor'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
        if 'EWUTil' in self.dict_data and 'ELEC' in self.dict_data:
            nstep = len(self.dict_data['ELEC'])
            for istep in range(nstep):
                self.dict_data['ELEC'][istep] = self.dict_data['ELEC'][istep] + self.dict_data['EWUTil'][istep]
                self.dict_text['ELEC'][istep] = "%f" % self.dict_data['ELEC'][istep]
            
    def test_diff(self, obj, tolerance): # compare energies
        # test MD steps
        is_failure = False
        dict_failure = dict.fromkeys(list(self.dict_data.keys()), False)
        nstep_failure = 0

        nstep = len(self.dict_data['Step'])
        for istep in range(nstep):
            for key in list(self.dict_data.keys()):
                d = abs(self.dict_data[key][istep] - obj.dict_data[key][istep])
                if d > tolerance:
                    is_failure = True
                    dict_failure[key] = True
                    nstep_failure = istep
            if is_failure:
                break

        if is_failure:
            self.is_passed = False
            print("Failure at step %d (tolerance = %4.2e)" % (self.dict_data['Step'][nstep_failure], tolerance))
            nstep_max = min([nstep_failure + 3, nstep])
            for istep in range(nstep_failure, nstep_max):
                print("Step %d" % (self.dict_data['Step'][istep]))

                sys.stdout.write("  ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % key.rjust(14))
                sys.stdout.write("\n")

                sys.stdout.write("< ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % self.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n")
                
                sys.stdout.write("> ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % obj.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n\n")
        else:
            self.is_passed = True
            print("Passed (tolerance = %4.2e)" % (tolerance))
            
    def test_diff_energies(self, obj, tolerance): # compare energies ignoring HFCKe and PRESSE terms
        # test MD steps
        is_failure = False
        dict_failure = dict.fromkeys(list(self.dict_data.keys()), False)
        nstep_failure = 0

        nstep = len(self.dict_data['Step'])
        for istep in range(nstep):
            for key in list(self.dict_data.keys()):
                if ((key == "Time") or (key == "TOTEner") or (key == "TOTKe") or (key == "ENERgy") or (key == "TEMPerature") or (key == "GRMS") or (key == "BONDs") or (key == "ANGLes") or (key == "UREY-b") or (key == "DIHEdrals") or (key == "IMPRopers") or (key == "CMAPs") or (key == "VDWaals") or (key == "ELEC")):
                    d = abs(self.dict_data[key][istep] - obj.dict_data[key][istep])
                    if d > tolerance:
                        is_failure = True
                        dict_failure[key] = True
                        nstep_failure = istep
            if is_failure:
                break

        if is_failure:
            self.is_passed = False
            print("Failure at step %d (tolerance = %4.2e)" % (self.dict_data['Step'][nstep_failure], tolerance))
            nstep_max = min([nstep_failure + 3, nstep])
            for istep in range(nstep_failure, nstep_max):
                print("Step %d" % (self.dict_data['Step'][istep]))

                sys.stdout.write("  ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % key.rjust(14))
                sys.stdout.write("\n")

                sys.stdout.write("< ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % self.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n")
                
                sys.stdout.write("> ")
                for key in list(self.dict_data.keys()):
                    if dict_failure[key]:
                        sys.stdout.write("%14s" % obj.dict_text[key][istep].rjust(14))
                sys.stdout.write("\n\n")
        else:
            self.is_passed = True
            print("Passed (tolerance = %4.2e)" % (tolerance))
            
    def calc_relative_error(self):
        print("not available.")

############### TEST ##########################################################
if __name__ == '__main__':

    # test of 13 digits charmm format generated by GENESIS
    text = '''[STEP4] Compute Single Point Energy for Molecules
 
Output_Energy> CHARMM_Style is used
 
DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        0      0.00000      0.00000      0.00000      0.00000      0.00000
DYNA PROP>          0.00000      0.00000      0.00000      0.00000      0.00000
DYNA INTERN>      585.70255   1228.72913    167.39536    921.88694    102.07776
DYNA EXTERN>     7826.04840 -82787.99150      0.00000      0.00000      0.00000
DYNA IMAGES>        0.00000      0.00000      0.00000      0.00000      0.00000
DYNA EWALD>         0.00000      0.00000      0.00000      0.00000      0.00000
DYNA PRESS>         0.00000      0.00000      0.00000      0.00000      0.00000
 ----------       ---------    ---------    ---------    ---------    ---------
 
[STEP5] Perform Molecular Dynamics Simulation
 
Initial_Velocity> Generate initial velocities
  iseed           =  113420321
  temperature     =    298.150
 
Stop_Trans_Rotation> Information about center of mass
  position        =     32.818221        32.586474        33.117938    
  velocity        =   -0.10577230E-02  -0.21867632E-02   0.12647128E-02
  angul_momentum  =    -2103.9577       -8227.0787       -709.00038    
  kinetic_ene     =    0.54193315    
 
Stop_Trans_Rotation> Translational motion was removed
 
Output_Energy> CHARMM_Style is used
 
DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        0      0.00000 -57645.31039  14291.91692 -71937.22731    297.30761
DYNA PROP>         13.96023      0.00000  14325.47320      0.00000  12961.05092
DYNA INTERN>      587.17584   1227.87650    156.01780    921.78016    101.45686
DYNA EXTERN>     7830.32767 -82761.86215      0.00000      0.00000      0.00000
DYNA PRESS>         0.00000  -8640.70061   2711.92949    252.53539 240990.21157
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        1      0.00100 -57642.77801  14292.39734 -71935.17535    297.31760
DYNA PROP>         13.95558      0.00000  14325.86224      0.00000  12692.20822
DYNA INTERN>      586.23026   1230.73940    156.30088    920.53865    103.33305
DYNA EXTERN>     7822.17684 -82754.49444      0.00000      0.00000      0.00000
DYNA PRESS>         0.00000  -8461.47215   2712.02065    303.64017 240990.21157
 ----------       ---------    ---------    ---------    ---------    ---------
'''

    ene = Charmm()

    # test parse
    ene.parse(text.split("\n"))
    assert len(ene.dict_data) == 26
    assert ene.dict_text['ELEC'][0].lstrip() == "-82761.86215"
    assert ((ene.dict_data['ELEC'][0] - 82761.86215) < 10**(-3))

    # test diff
    ene2 = copy.deepcopy(ene)
    ene2.dict_text['ELEC'][0] = "-82762.92832"
    ene2.dict_data['ELEC'][0] = -82762.92832
    ene2.dict_text['DIHEdrals'][0] = "921.79238"
    ene2.dict_data['DIHEdrals'][0] = 921.79238
    ene.test_diff(ene2, 0.0001)

    # test of 14 digits charmm format generated by GENESIS
    text = '''[STEP4] Compute Single Point Energy for Molecules
 
Output_Energy> CHARMM_Style is used
 
DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA CROSS:           CMAPs                                                    
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
DYNA EWALD:          EWKSum       EWSElf       EWEXcl       EWQCor       EWUTil
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        0       0.00000       0.00000       0.00000       0.00000       0.00000
DYNA PROP>           0.00000       0.00000       0.00000       0.00000       0.00000
DYNA INTERN>     50626.48594   41901.36860     718.14958   13934.58410     369.13107
DYNA CROSS>      -4835.31122
DYNA EXTERN>    180164.31137-1681535.91075       0.00000       0.00000       0.00000
DYNA IMAGES>         0.00000       0.00000       0.00000       0.00000       0.00000
DYNA EWALD>          0.00000       0.00000       0.00000       0.00000       0.00000
DYNA PRESS>          0.00000       0.00000       0.00000       0.00000       0.00000
 ----------       ---------    ---------    ---------    ---------    ---------
 
[STEP5] Perform Molecular Dynamics Simulation
 
Initial_Velocity> Generate initial velocities
  iseed           =     314159
  temperature     =    303.000
 
Stop_Trans_Rotation> Information about center of mass
  position        =     90.975571        73.945738        69.101611    
  velocity        =   -0.10096501E-02  -0.68857264E-03  -0.41584577E-03
  angul_momentum  =     33445.157       -74978.687       -51807.025    
  kinetic_ene     =     1.6889621    
 
Stop_Trans_Rotation> Translational motion was removed
Stop_Trans_Rotation> Rotational motion was removed
 
Output_Energy> CHARMM_Style is used
 
DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA CROSS:           CMAPs                                                    
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        0       0.00000-1102973.88152  295683.30978-1398657.19130     303.08506
DYNA PROP>           0.68936       0.00000  295692.24304       0.00000   58224.74527
DYNA INTERN>     50626.48594   41901.36860     718.14958   13934.58410     369.13107
DYNA CROSS>      -4835.31122
DYNA EXTERN>    180164.31137-1681535.91075       0.00000       0.00000       0.00000
DYNA PRESS>          0.00000  -38816.49684       0.00000    3497.82827 3104392.80981
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>       10       0.01000-1102171.72790  199241.40716-1301413.13506     204.22896
DYNA PROP>           6.73354       0.00000  200056.67939       0.00000 -231313.62819
DYNA INTERN>     50275.97200   52933.75308    2489.77467   15630.61920    1476.24534
DYNA CROSS>      -4679.64595
DYNA EXTERN>    191984.83464-1611524.68804       0.00000       0.00000       0.00000
DYNA PRESS>          0.00000  154209.08546       0.00000    6342.19102 3104392.80981
 ----------       ---------    ---------    ---------    ---------    ---------
Update_Pairlist_Pbc> Search Nonbonded Interactions 
within PairListDist
  num_nb15        =  124855668
 
DYNA>       20       0.02000-1101503.89915  177696.55379-1279200.45294     182.14478
DYNA PROP>           8.59759       0.00000  179193.42738       0.00000 -320688.59931
DYNA INTERN>     59219.07638   59489.86932    2503.37341   15329.32038     802.20635
DYNA CROSS>      -4618.48371
DYNA EXTERN>    215022.87035-1626948.68543       0.00000       0.00000       0.00000
DYNA PRESS>          0.00000  213792.39954       0.00000    7341.34682 3104392.80981
 ----------       ---------    ---------    ---------    ---------    ---------
'''

    ene = Charmm()

    # test parse
    ene.parse(text.split("\n"))
    assert len(ene.dict_data) == 27
    assert ene.dict_text['ELEC'][0].lstrip() == "-1681535.91075"

    # test diff
    ene2 = copy.deepcopy(ene)
    ene2.dict_text['ELEC'][0] = "-1681535.93333"
    ene2.dict_data['ELEC'][0] = -1681535.93333
    ene.test_diff(ene2, 0.0001)

    # test of original charmm format generated by CHARMM
    text = ''' CHARMM>    !dynamics LEAP VERLET nstep 1 timestp 0.001 -
 CHARMM>    dynamics VV2 VERLET nstep 1 timestp 0.001 -
 CHARMM>             firstt 0.0 tstruc 0.0 iseed 314159 -
 CHARMM>             inbfrq 1 nprint 1
  IUNREA = -1         IUNWRI = -1          IUNOS = -1
  IUNCRD = -1         IUNVEL = -1          KUNIT = -1

 NONBOND OPTION FLAGS: 
     ELEC     VDW      ATOMs    CDIElec  SHIFt    VATOm    VSWItch 
     BYGRoup  NOEXtnd  NOEWald 
 CUTNB  = 12.000 CTEXNB =999.000 CTONNB =  8.000 CTOFNB = 10.000
 WMIN   =  1.500 WRNMXD =  0.500 E14FAC =  1.000 EPS    =  1.000
 NBXMOD =      5
 There are  6401815 atom  pairs and    34709 atom  exclusions.
 There are        0 group pairs and     2149 group exclusions.
   NSTEP =        1    NSAVC =       10    NSAVV =       10
  ISCALE =        0   ISCVEL =        0   IASORS =        0
  IASVEL =        1   ICHECW =        1   NTRFRQ =        0
  IHTFRQ =        0   IEQFRQ =        0   NPRINT =        1
  INBFRQ =        1   IHBFRQ =        0   IPRFRQ =        1
  ILBFRQ =       50   IMGFRQ =       50    ISEED =              314159
  ISVFRQ =        0   NCYCLE =        1    NSNOS =       10
  FIRSTT =     0.000  TEMINC =     5.000  TSTRUC =     0.000
  FINALT =   298.150  TWINDH =    10.000  TWINDL =   -10.000

  TIME STEP =  2.04548E-02 AKMA       1.00000E-03 PS


           SHAKE TOLERANCE =     0.10000E-07
 NUMBER OF DEGREES OF FREEDOM =  48387

          SEED FOR RANDOM NUMBER GENERATOR IS       314159
          GAUSSIAN OPTION                  IS            1
          VELOCITIES ASSIGNED AT TEMPERATURE =      0.0000

     DETAILS ABOUT CENTRE OF MASS
     POSITION          :  -3.67788895E-02  -7.52627313E-03  -9.06203823E-03
     VELOCITY          :    0.0000000        0.0000000        0.0000000    
     ANGULAR MOMENTUM  :    0.0000000        0.0000000        0.0000000    
     KINETIC ENERGY    :    0.0000000    
DYNA DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
DYNA PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
DYNA INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
DYNA EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
DYNA IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
DYNA PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        0      0.00000 -75725.03876      0.00000 -75725.03876      0.00000
DYNA PROP>         14.57730      0.00000      0.00000      0.00000   9793.81757
DYNA INTERN>      587.17584   1227.87650    156.01780    921.78016    101.45686
DYNA EXTERN>     6665.06115 -74035.30372      0.00000      0.00000      0.00000
DYNA IMAGES>      998.43105 -12347.53442      0.00000      0.00000      0.00000
DYNA PRESS>      -766.16674  -5763.04497    217.99525  -1639.74283 240990.21157
 ----------       ---------    ---------    ---------    ---------    ---------
DYNA>        1      0.00100 -75726.15470    130.81846 -75856.97316      2.72101
DYNA PROP>         14.52245      0.00000      0.00000      0.00000   9790.85358
DYNA INTERN>      570.25760   1198.54478    151.55814    919.37130     99.50939
DYNA EXTERN>     6655.52114 -74091.33853      0.00000      0.00000      0.00000
DYNA IMAGES>      997.13467 -12357.53166      0.00000      0.00000      0.00000
DYNA PRESS>      -630.42946  -5896.80626    179.37431  -1652.98728 240990.21157
 ----------       ---------    ---------    ---------    ---------    ---------

 * * * AVERAGES FOR THE LAST        1 STEPS
AVER DYN: Step         Time      TOTEner        TOTKe       ENERgy  TEMPerature
AVER PROP:             GRMS      HFCTote        HFCKe       EHFCor        VIRKe
AVER INTERN:          BONDs       ANGLes       UREY-b    DIHEdrals    IMPRopers
AVER EXTERN:        VDWaals         ELEC       HBONds          ASP         USER
AVER IMAGES:        IMNBvdw       IMELec       IMHBnd       RXNField    EXTElec
AVER PRESS:            VIRE         VIRI       PRESSE       PRESSI       VOLUme
 ----------       ---------    ---------    ---------    ---------    ---------
AVER>        1      0.00100 -75726.15470    130.81846 -75856.97316      2.72101
AVER PROP>         14.52245      0.00000      0.00000      0.00000   9790.85358
AVER INTERN>      570.25760   1198.54478    151.55814    919.37130     99.50939
AVER EXTERN>     6655.52114 -74091.33853      0.00000      0.00000      0.00000
AVER IMAGES>      997.13467 -12357.53166      0.00000      0.00000      0.00000
AVER PRESS>      -630.42946  -5896.80626    179.37431  -1652.98728 240990.21157
 ----------       ---------    ---------    ---------    ---------    ---------

 * * * RMS FLUCTUATIONS FOR         1 STEPS
FLUC>        1          0.0 5734450506.1      17113.5 5754280377.7          7.4
FLUC PROP>          210.901        0.000        0.000        0.000 95860813.856
FLUC INTERN>    325193.7255 1436509.5936   22969.8713  845243.5857    9902.1182
FLUC EXTERN>     44295961.7 5489526445.2          0.0          0.0          0.0
FLUC IMAGES>      994277.55 152708588.73         0.00         0.00         0.00
FLUC PRESS>         397441.    34772324.       32175.     2732367. 58076282071.
 ----------       ---------    ---------    ---------    ---------    ---------

     DETAILS ABOUT CENTRE OF MASS
     POSITION          :  -3.67788895E-02  -7.52627313E-03  -9.06203823E-03
     VELOCITY          :  -1.36816176E-15  -3.16622913E-16  -2.71919912E-16
     ANGULAR MOMENTUM  :    4.6459702        138.80512       -45.385133    
     KINETIC ENERGY    :   1.47839360E-25
  
 CHARMM>    stop
$$$$$$  New timer profile $$$$$
   Shake Setup                     0.01 Other:            0.00
   First List                      2.94 Other:            0.00
   Shake time                      0.02 Other:            0.00
   Dynamics total                  0.01 Other:            0.00
         Electrostatic & VDW             0.63 Other:            0.00
      Nonbond force                   0.63 Other:            0.00
         Bond energy                     0.00 Other:            0.00
         Angle energy                    0.00 Other:            0.00
         Dihedral energy                 0.00 Other:            0.00
         Restraints energy               0.00 Other:            0.00
      INTRNL energy                   0.01 Other:            0.00
   Energy time                     0.64 Other:            0.00
 Total time                      8.61 Other:            4.99

                    NORMAL TERMINATION BY NORMAL STOP
                    MAXIMUM STACK SPACE USED IS 9478830
                    STACK CURRENTLY IN USE IS         0
                    MOST SEVERE WARNING WAS AT LEVEL  1
                    HEAP PRINTOUT-  HEAP SIZE         10240000
                    SPACE CURRENTLY IN USE IS            13824
                    MAXIMUM SPACE USED IS             29464250
                    FREE LIST
  PRINHP> ADDRESS:               1 LENGTH:        10225676 NEXT:        10239501
  PRINHP> ADDRESS:        10239501 LENGTH:             500 NEXT:       199443967
  PRINHP> ADDRESS:       199443967 LENGTH:          671744 NEXT:       201278975
  PRINHP> ADDRESS:       201278975 LENGTH:        32309248 NEXT:               0

                    $$$$$ JOB ACCOUNTING INFORMATION $$$$$
                     ELAPSED TIME:     8.62  SECONDS 
                         CPU TIME:     8.58  SECONDS 
'''

    ene = Charmm()

    # test parse
    ene.parse(text.split("\n"))
    assert len(ene.dict_data) == 26
    print(ene.dict_text['ELEC'][0])
    assert ene.dict_text['ELEC'][0] == "-86382.838140"

    # test diff
    ene2 = copy.deepcopy(ene)
    ene2.dict_text['ELEC'][0] = "-86382.923020"
    ene2.dict_data['ELEC'][0] = -86382.923020
    ene.test_diff(ene2, 0.0001)

