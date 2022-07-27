#!/usr/bin/env python
# -*- encoding: utf-8 -*-

## @package fortdep
#  @brief   fortran dependency inspector
#  @author  Motoshi Kamiya
#
#  (c) Copyright 2016-2017 RIKEN. All rights reserved.

# Simple fortran dependency checker for GENESIS
#
# usage: fortdep.py [list of files (*.fpp)]

from __future__ import print_function
import os
import re
import sys
import getopt
import glob

class FortranFile:

  # class variables
  re_fcomment = re.compile( "(^[^!]*)!(.*$)" )
  re_module   = re.compile( "^(.*;+|\s*)\s*module\s*([^\s,]*)\s*", re.I )
  re_use      = re.compile( "^(.*;+\s*|\s*)use\s+([^\s,]*)\s*", re.I )
  mod_ext     = ".mod"

  def __init__( self, fname = "", ext = ".o" ):
    self.setFilename( fname, ext )
    self.filename = fname
    self.modules  = []
    self.depmods  = []

  def setFilename( self, fname, ext = ".o" ):
    if len(fname) == 0:
      self.filename = ""
      self.objname  = ""
      return
    self.filename = fname
    self.objname = re.sub( r'\.[a-zA-Z0-9]+$', ext, self.filename )

  def parse( self ):
    myf = open( self.filename, 'r' )
    for line in myf:
      ma = FortranFile.re_fcomment.search( line )
      if ( ma ): # replace string
        line = ma.group(1)
      mb = FortranFile.re_module.search( line )
      if ( mb ):
        modname = mb.group(2).lower()
        self.appendModuleNames( modname )
      ma = FortranFile.re_use.search( line )
      if ( ma ):
        modname = ma.group(2).lower()
        self.appendDependMods( modname )
    myf.close()
    self.makeUniqList()

  def appendModuleNames( self, modname ):
    self.modules.append( modname + FortranFile.mod_ext )

  def appendDependMods( self, modname ):
    self.depmods.append( modname + FortranFile.mod_ext )

  def makeUniqList( self ):
    self.modules = list(set(self.modules))
    self.depmods = list(set(self.depmods))

  def getMyModuleFilenames( self ):
    return self.modules

  def getDepModuleFilenames( self ):
    return self.depmods

  def showModuleList( self, output ):
    for m in self.modules:
      print( m, file=output )

  def showDepModList( self, output ):
    for m in self.depmods:
      print( m, file=output )

  def recipe( self, mods_avail = [], static_deps = "" ):
    depmods = []
    for m in self.depmods:
      if m.lower() in mods_avail or not mods_avail:
        depmods.append(m)
    ret = ""
    ret += self.objname + ": " + self.filename + " " + " ".join(depmods) + " " + static_deps
    if len(self.modules) > 0:
      ret += "\n"
      ret += " ".join(self.modules) + ": " + self.filename + " " + self.objname
    return ret

def usage( ret = 0 ):
  print( "fortdep: fortran dependency inspector", file = sys.stderr )
  print( "usage:   fortdep [options] files > [output]", file = sys.stderr )
  print( "options:", file = sys.stderr )
  print( "        -s [static dependencies for obj]", file = sys.stderr )

  sys.exit(ret)

if __name__ == "__main__":
  static_deps = ""

  # for future extension
  try:
    opts, args = getopt.getopt( sys.argv[1:], "hs:e:f:" )
  except getopt.GetoptError:
    print( "Error, failed to parse options", file = sys.stderr )

  if len(args) == 0:
    usage()

  # parse opts
  for o,a in opts:
    if o in ( "-s" ):
      static_deps = a
    elif o in ( "-h" ):
      usage()

  files = []
  mods_in_this_dir = []
  # build a list of modules
  for f in args:
    ff = FortranFile( f )
    ff.parse()
    files.append( ff )
    mods = ff.getMyModuleFilenames()
    for m in mods:
      mods_in_this_dir.append( m )

  # unique
  mods_in_this_dir = list(set(mods_in_this_dir))

  ## for debug
  #for ff in files:
  #  ff.showModuleList( sys.stderr )
  #  ff.showDepModList( sys.stderr )

  for ff in files:
    print( ff.recipe( mods_in_this_dir, static_deps ) )
