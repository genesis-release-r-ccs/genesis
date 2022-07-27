#!/usr/bin/env python

#
# A python script for GENESIS analysis tool regression tests
#
# (c) Copyright 2014 RIKEN. All rights reserved.
#

# Usage:
#
#    $ ./test_analysis.py /home/user/genesis/bin/
#
# With the default settings, this script scans the directory where
# it is located for folders starting with "test_" and generates
# a separate test class for each folder.
# The analysis program to run is determined from the folder name
# after removing the "test_" prefix.
# Each test class then searches for folders that contain no subdirectories
# and generates a separate test for each one of them.
# By default, a test folder must have an "int" input script
# and a "ref" file with the expected analysis output.
# An "out" output file must be produced after the execution
# of the analysis program inside the test folder,
# which is then compared with the "ref" file.
# Divergence tolerance, analysis program execution parameters
# and environmental variables can be set in "config.ini" files
# located inside the test folder or any of its parent folders.
# When same options are set in several folders,
# values set in children folders overwrite those set in parent folders.
# The "config.ini" format is identical to that of "defaults.ini".
# By default, the script will only search the test folder for
# analysis program executables, so directories containing
# analysis programs should be provided via command line arguments
# or the program_folder option in the configuration files.

import unittest
import sys
import os
import os.path
import subprocess

# python3 compatibility
try:
    from ConfigParser import SafeConfigParser as ConfigParser
except ImportError:
    from configparser import ConfigParser


# subclass ConfigParser to preserve case in option names
class CaseConfigParser(ConfigParser):

    def optionxform(self, option):
        return option

# a metaclass for TestAnalysisProgramBase that generates
# separate test_ methods for each bottommost directory under self.topdir
class TestAnalysisProgramMeta(type):

    def __new__(mcs, name, bases, dict):

        def generate_folder_test(folder):

            def test_folder(self):

                self.run_test_in_folder(folder)

            return test_folder

        # self.topdir will be defined dynamically

        if "topdir" in dict:

            topdir = dict["topdir"]

            # generate a test for each bottommost folder

            for path, folders, files in os.walk(topdir):

                if len(folders) == 0:

                    module_name = "test_folder_" + path

                    dict[module_name] = generate_folder_test(path)


        return type.__new__(mcs, name, bases, dict)

# python3 metaclass compatibility
try:
    exec("""class TestAnalysisProgramPreBase(
    metaclass=TestAnalysisProgramMeta): pass""")
except SyntaxError:
    class TestAnalysisProgramPreBase(object):
        __metaclass__ = TestAnalysisProgramMeta

class TestAnalysisProgramBase(TestAnalysisProgramPreBase):

    # folders where tests will look for analysis program executables
    program_folders = []

    # search for executables inside program_folders and current folder
    @classmethod
    def find_executable(cls, folder, config):

        section_settings = "settings"

        exe_name = config.get(section_settings, "program_name")

        search_paths = [folder] + cls.program_folders

        program_folder = config.get(section_settings, "program_folder")
        if program_folder:
            search_path.append(program_folder)

        for search_path in search_paths:

            search_path = os.path.abspath(search_path)

            path = os.path.join(search_path, exe_name)
            if os.path.isfile(path) and os.access(path, os.X_OK):
                return path

        raise Exception("Executable '%s' not found. " % (exe_name,)
              + "Did you set program_folder correctly by passing "
              + "folder with analysis programs to "
              + "%s or setting it in %s?" % (sys.argv[0], cls.defaults_ini))

    # compare numeric values in reference and output files
    def check_output_file(self, folder, config):

        section_name = "settings"
        ref_name = config.get(section_name, "reference_name")
        out_name = config.get(section_name, "output_name")

        ref_path = os.path.join(folder, ref_name)
        out_path = os.path.join(folder, out_name)

        ref_fobj = open(ref_path, "r")
        out_fobj = open(out_path, "r")

        ref_lines = ref_fobj.readlines()
        out_lines = out_fobj.readlines()

        ref_fobj.close()
        out_fobj.close()

        self.assertEqual(len(ref_lines), len(out_lines),
                "Output file '%s' " % (out_path,)
              + "has different number of lines than reference file")

        tolerance = config.get(section_name, "tolerance").split()
        tolerance = tuple(map(float, tolerance))
        minimal_tolerance = min(tolerance)

        for irow, (rline, oline) in enumerate(zip(ref_lines, out_lines)):
            rdata = rline.split()
            odata = oline.split()

            self.assertEqual(len(rdata), len(odata),
                    "Different number of columns in "
                  + "'%s' at row %d" % ( out_path, irow))

            for icol, (rstr, ostr) in enumerate(zip(rdata, odata)):

                try:
                    rval = float(rstr)
                except ValueError:
                    rval = None

                try:
                    oval = float(ostr)
                except ValueError:
                    oval = None

                if oval is None:
                    self.assertEqual(rval, oval,
                            "Could not convert to float at "
                          + "row %d, column %d" % (irow, icol)
                          + " in file '%s'" % (out_path,))
                    continue

                if rval is None:
                    continue

                try:
                    error = tolerance[icol]
                except IndexError:
                    error = minimal_tolerance

                diff = abs(oval - rval)
                self.assertTrue(diff <= error,
                        "Divergence larger than tolerance "
                      + "(|%g - %g| = %g > %g) " % (oval, rval, diff, error)
                      + "at row %d, column %d " % (irow, icol)
                      + "in file '%s'" % (out_path,))


    @classmethod
    def get_configuration(cls, folder):

        # cls.defaults_ini will be appended dynamically
        defaults = CaseConfigParser()
        defaults.read(cls.defaults_ini)

        config_name = defaults.get("settings", "config_name")

        # store all possible configuration file paths in config_inis
        config_inis = []
        path = folder
        while True:
            config_inis.insert(0, os.path.join(path, config_name))
            if os.path.samefile(path, cls.topdir): break
            path = os.path.dirname(path)
        config_inis.insert(0, cls.defaults_ini)

        config = CaseConfigParser()
        config.read(config_inis)

        # set program name if blank
        section_name = "settings"
        option_name = "program_name"

        test_prefix = config.get(section_name, "test_folder_prefix")

        # determine program name from topdir if not explicitly set
        if not config.get(section_name, option_name):
            program_name = os.path.basename(cls.topdir)
            if program_name.startswith(test_prefix):
                program_name = program_name[len(test_prefix):]
            config.set(section_name, option_name, program_name)

        return config

    # check if input and reference files are inside folder
    @staticmethod
    def check_test_folder(folder, config):

        section_name = "settings"
        option_names = ["input_name", "reference_name"]

        for option_name in option_names:
            file_name = config.get(section_name, option_name)
            path = os.path.join(folder, file_name)
            if not os.path.isfile(path):
                raise Exception(
                        "File '%s' does not exist in %s" % (file_name, folder))

    # remove output files from test folder, because most GENESIS analysis
    # programs can't overwrite them
    @staticmethod
    def prepare_test_folder(folder, config):

        out_name = config.get("settings", "output_name")
        out_path = os.path.join(folder, out_name)

        if os.path.exists(out_path): os.remove(out_path)


    def execute_analysis_program_in_folder(self, folder, config):

        exe = self.find_executable(folder, config)

        # set program environment according to config
        env = os.environ.copy()
        section_env = "environ"
        for option in config.options(section_env):

            # delete option if value is an empty string,
            # set to vale otherwise
            value = config.get(section_env, option)
            if value:
                env[option] = config.get(section_env, option)
            else:
                if option in env: del env[option]

        section_settings = "settings"
        inp_name = config.get(section_settings, "input_name")
        log_name = config.get(section_settings, "stdout_name")
        err_name = config.get(section_settings, "stderr_name")

        log_path = os.path.join(folder, log_name)
        err_path = os.path.join(folder, err_name)

        log_fobj = open(log_path, "w")
        err_fobj = open(err_path, "w")

        args = [exe, inp_name]

        mpi_name = config.get(section_settings, "mpi_name")

        # use mpi execution if mpi_name is not an empty string
        if mpi_name:

            mpi_num = config.get(section_settings, "mpi_process_number")

            args[0:0] = [mpi_name, "-n", mpi_num]

        proc = subprocess.Popen(
                args,
                stdout=log_fobj, stderr=err_fobj,
                cwd=folder,
                env=env)

        ret = proc.wait()

        log_fobj.close()
        err_fobj.close()

        self.assertEqual(ret, 0, """Analysis program did not terminate normally.
Execution arguments: %s\n""" % (" ".join(args),))


    # this method will be called by all dynamically generated test_ methods
    def run_test_in_folder(self, folder):

        config = self.get_configuration(folder)

        self.check_test_folder(folder, config)

        self.prepare_test_folder(folder, config)

        self.execute_analysis_program_in_folder(folder, config)

        self.check_output_file(folder, config)

# generates separate test class for each test_ folder
def generate_test_classes():

    # use the folder where test.py is located
    work_dir = os.path.dirname(os.path.abspath(os.path.realpath(__file__)))

    # default values are stored in default.ini
    defaults_ini = os.path.join(work_dir, "defaults.ini")

    defaults = CaseConfigParser()
    defaults.read(defaults_ini)

    test_prefix = defaults.get("settings", "test_folder_prefix")

    # create separate class for each folder
    path, folders, files = next(os.walk(work_dir))

    for folder in folders:
        if not folder.startswith(test_prefix): continue
        test_name = folder[len(test_prefix):]
        class_name = "Test_" + test_name
        globals()[class_name] = type(class_name,
                (TestAnalysisProgramBase, unittest.TestCase),
                dict(
                    topdir=os.path.join(path, folder),
                    defaults_ini=defaults_ini,
                    ))

generate_test_classes()
# delete just in case, so it can never be called again
del generate_test_classes

def main(args=None):

    import optparse

    usage = """Usage: %prog [-h] [program_folder ...] [-- unittest_arguments ...]
Unittest help: %prog -- -h

Positional arguments:
  program_folder      folder where analysis programs are stored"""

    parser = optparse.OptionParser(usage)

    # arguments before "--" are for optparse
    # arguments after "--" are for unittest
    if args is None:
        try:
            idx = sys.argv.index("--")
        except ValueError:
            idx = len(sys.argv)

        # arguments for optparse
        args = sys.argv[1:idx]

        # arguments for unittest
        sys.argv = [sys.argv[0]] + sys.argv[idx+1:]

        # make varbose by default
        if len(sys.argv) == 1: sys.argv.append("-v")

    (options, args) = parser.parse_args(args)

    for folder in args:
        TestAnalysisProgramBase.program_folders.append(
                os.path.abspath(os.path.realpath(folder)))

    unittest.main()

if __name__ == "__main__":
    main()
