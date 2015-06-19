# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 14:33:52 2012

@author: jim
"""

import os
import subprocess
import sys

def check_create_dir( strDir ):

	if not os.path.exists( strDir ):
		os.makedirs( strDir )
	return strDir

def check_file(strPath):
	try:
		with open(strPath):
			pass
	except IOError:
		print( strPath, "does not exist. Please check the path to this file." )
		return

def CheckDependency(strCmd,strArg,strIntendedProgram):

    if (os.path.isfile(strCmd)==False):
        raise IOError("\nShortBRED was unable to find " +  strIntendedProgram + " at the path *" + strCmd +  "*\nPlease check that the program is installed and the path is correct.")
        #print "\nShortBRED was unable to load " +  strIntendedProgram + "at the path " + strCmd +  "\nPlease check to make sure the program is installed and the path is correct."
        sys.exit(1)

    if (strArg!=""):
        pCmd = subprocess.Popen([strCmd,strArg],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    else:
        pCmd = subprocess.Popen([strCmd],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    outCmd= pCmd.communicate()[0]
    if strIntendedProgram == "cdhit":
        sys.stderr.write("Path for cdhit appears to be fine. This program returns an error [exit code=1] when tested and working properly, so ShortBRED does not check it.\n")
    else:
        if( pCmd.returncode==0):
            sys.stderr.write("Tested " + strIntendedProgram  + ". Appears to be working.\n")
        else:
            sys.stderr.write("Tested " + strIntendedProgram  + " returned a nonzero exit code (typically indicates failure). Please check to ensure the program is working. Will continue running.\n")