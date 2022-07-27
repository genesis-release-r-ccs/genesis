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
class Fep(object):
    def __init__(self):  # initialization
        # public attributes
        self.dict_text = {}
        self.dict_data = {}
        self.dict_error = {}
        self.is_passed = False

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
        title = []
        text = []
        data = []
        data_each = []
        is_md = False
        for line in text_list:
            if not line.startswith("#"):
                data_each = line.rstrip('\n').split()
                if len(data_each) > 0:
                    text.append(data_each)
                    tmp = []
                    for i in range(len(data_each)):
                        if data_each[i][-1].isdigit(): # number ?
                            tmp = tmp + [(float(data_each[i]))]
                        else: # Not a Number !
                            tmp = tmp + [(float(999999))] 
                    data.append(tmp)
            else:
                result = re.search("STEP", line)
                if result is not None:
                    line_split = line.rstrip('\n').strip('#').split()
                    title = title + line_split

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

    def test_diff(self, obj, tolerance): # compare energies
        # test MD steps
        keys = set(self.dict_data.keys()) & set(obj.dict_data.keys())
        is_failure = False
        dict_failure = dict.fromkeys(keys, False)
        nstep_failure = 0

        nstep = len(self.dict_data['STEP'])
        for istep in range(nstep):
            for key in keys:
                d = abs(self.dict_data[key][istep] - obj.dict_data[key][istep])
                if abs(self.dict_data[key][istep]) < 1e4: #min log is 1e-4
                    ratio=d
                else:
                    ebase=max(abs(self.dict_data[key][istep]),1.0)
                    ratio = d/ebase
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
            print("Failure at step %d (tolerance = %4.2e(ene))" % (self.dict_data['STEP'][nstep_failure], tolerance))
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
            print("Passed (tolerance = %4.2e(ene))" % (tolerance))

############### TEST ##########################################################
if __name__ == '__main__':
    text = '''# FEP window index        1
#     STEP     Total_E_ref     Delta_E_rev     Delta_E_fwd
        55     -30601.3858         -0.0000        -11.9887
        60     -30601.5122          0.0000        -12.4823
        65     -30593.9519          0.0000        -12.5749
        70     -30602.1777          0.0000        -11.6161
        75     -30600.3690         -0.0000        -11.8543
        80     -30593.0163          0.0000        -12.7094
        85     -30594.2356          0.0000        -12.3298
        90     -30597.9725         -0.0000        -11.4078
        95     -30598.5747          0.0000        -12.0962
       100     -30591.5799          0.0000        -13.1131
# FEP window index        2
       155     -30593.4023         14.4644         -0.0000
       160     -30592.9836         13.8435         -0.0000
       165     -30593.3225         13.7767         -0.0000
       170     -30592.0947         15.1606         -0.0000
       175     -30586.4292         15.1691          0.0000
       180     -30577.7778         14.0791          0.0000
       185     -30580.2632         13.5071         -0.0000
       190     -30574.9835         13.9099         -0.0000
       195     -30580.5005         13.1632          0.0000
       200     -30582.4343         12.9324         -0.0000
'''
    ene = Fep()

    # test parse
    ene.parse(text.split("\n"))
#     assert len(ene.dict_data) == 25
#     assert ene.dict_text['ELECT'][0] == "-39052.9065"
#     assert ((ene.dict_data['ELECT'][0] + 39052.9065) < 10**(-3))

    # test diff
    ene2 = copy.deepcopy(ene)
    ene2.dict_text['Delta_E_rev'][1] = "-0.3"
    ene2.dict_data['Delta_E_rev'][1] = -0.3
    ene.test_diff(ene2, 1.0e-8)

    ene.delete_last()
