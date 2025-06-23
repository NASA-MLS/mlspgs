#!/usr/bin/env python3
from __future__ import print_function
import os, sys, re
import optparse as op
import subprocess as sp


#import glob

# this routine depends on teh IDL runtime saveset
# checknrtchunkforbreaks.sav. To create this runtime, do the following.
#
# Create a directory where you want to build the runtime saveset.
# copy over te file go to
# ~whdaffer/devel/mls/nrt/idl/make-checkchunkforbreaks-rt-saveset.pro
# to that directory.
#
# Make sure you have ~whdaffer/devel in your IDL path, with the '+'
# in front, so that IDL knows to traverse the tree while running the
# software to create the runtime.
#
# CD into that directory and do the following (the '%' is the sign for
# the shell prompt, so I want you to do this at the shell prompt.
#
# % idl make-checkchunkforbreaks-rt-saveset.pro
# 
# When that finishes, the runtime will be in the directory 
#
# Move that file to whatever you need to run the routine. Be aware
# that the working directory of an IDL runtime is the directory in
# which the runtime resides, *NOT* the working directory of the
# process that calls the runtime. If you need to run in a directory
# different from the one where the runtime lives, you can either copy
# it over to your working directory or set a symbolic link. 
#
# This routine makes the assumption that the runtime is in the CWD
# wheter this script is invoked.
# 
#

if __name__ == '__main__':

  usage='''%prog [options] file


This script is a simple wrapper for the idl routine
checknrtchunkforbreaks.pro. This routine takes the L2GP file for one
NRT chunk and checks it for `breaks'. If an error occurs, the IDL
exits with a status of -1, if the file has no breaks, IDL exits with a
status of 0 and if the file has 1 or more breaks, the status is 1

'''

  p=op.OptionParser(usage=usage)
  p.add_option('-v','--verbose',
               help='''emit more messages''',
               action='store_true',
               dest='verbose',
               default=False)

  p.add_option('-n','--nocatch',
               help='''skip catchblocks in IDL code''',
               action='store_true',
               dest='nocatch',
               default=False)


  print(f"python version: {sys.version}")
  opts,args=p.parse_args()
  nocatch=opts.nocatch
  verbose=opts.verbose

  if len(args) == 0:
    print(usage)
    print("Need file!")
    sys.exit(-1)


  # file=args[0]
  # if not os.path.isfile(file):
  #   print file + " doesn't exist!"
  #   sys.exit(-1)

  filename=args[0]
  spargs=['idl',
        '-rt=./checknrtgranuleforbreaks.sav',
        '-args',
        filename]

  if nocatch:
    spargs.append('nocatch=1')
  if verbose:
    spargs.append('verbose=1')
    print("Calling IDL with the following arguments")
    print(" ".join(spargs))
  
  pp=sp.Popen(spargs,stdout=sp.PIPE,stderr=sp.PIPE)
  (out,err)=pp.communicate()
  retcode=pp.returncode
  print("=================== STDOUT =============== ")
  print(str(out,'utf-8'))
  print( "================== STDERR =============== ")
  print( str(err,'utf-8'))
  print( "================== Done =============== ")

  if retcode == -1:
    print("Error in the IDL runtime!")
    sys.exit(1)
  elif retcode == 1:
    if verbose: print( "File has break!")
  elif retcode == 0:
    if verbose: 
      print( "File is okay!")
      print("retcode = {retcode:d}")
  else:
    print( "Some other error occurred")
    print(f"Error code = {retcode}")
    
  sys.exit(retcode)

          

  
#
# $Id$
#
# $Name$
#
# Modifications:
# $Log$
# Revision 1.8  2022/09/15 00:11:17  whdaffer
# Print STDOUT/STDERR regardless of return status from process running
# the IDL runtime. Need more info, not less. Deal with bytecode variables.
#
# Revision 1.7  2021/05/17 18:08:37  whdaffer
# Added from __future__ import print_function
#
# Revision 1.6  2021/05/05 14:28:16  whdaffer
# python 3 upgrade
#
# Revision 1.5  2017/06/26 22:31:30  whdaffer
# Took out dependency on verbose when the return status from the IDL is
# other than [-1, 0, 2]
#
# Revision 1.4  2017/06/26 21:12:52  whdaffer
# Added some more reportage to the `some other error occurred' message
#
# Revision 1.3  2017/05/18 23:06:29  whdaffer
# Added another branch to the test of return code from subprocess to
# handle any other error
#
# Revision 1.2  2017/02/23 19:34:24  whdaffer
# Changed idl -rt command line
#
# Revision 1.1  2017/02/23 19:33:49  whdaffer
# Renamed checknrtchunkforbreaks.py to checknrtgranuleforbreaks.py
#
# Revision 1.1  2017/02/16 22:53:41  whdaffer
# Initial revision
#
#                
# Copyright <year>, by the California Institute of
# Technology. ALL RIGHTS RESERVED. United States Government
# Sponsorship acknowledged. Any commercial use must be
# negotiated with the Office of Technology Transfer at the
# California Institute of Technology. 
#
# This software may be subject to U.S. export control
# laws. By accepting this software, the user agrees to
# comply with all applicable U.S. export laws and
# regulations. User has the responsibility to obtain export
# licenses, or other export authority as may be required
# before exporting such information to foreign countries or
# providing access to foreign persons. 
#
# Last Modified : Thu Apr 12 10:41:13 2012
