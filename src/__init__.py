# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 14:33:52 2012

@author: jim
"""

import os

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
