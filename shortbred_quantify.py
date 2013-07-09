#!/usr/bin/env python
#####################################################################################
#Copyright (C) <2013> Jim Kaminski and the Huttenhower Lab
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of ShortBRED (Short, Better REad Database)
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Jim Kaminski, jjk451@mail.harvard.edu).
#####################################################################################
import sys
import argparse
import subprocess
import csv
import re
import os
import datetime
import shutil
import tarfile
import gzip
import time
import math

import src
import src.quantify_functions
sq = src.quantify_functions

import numpy

import Bio
from Bio.Seq import Seq
from Bio import SeqIO


################################################################################
# Constants

c_strUSEARCH	= "usearch"

c_iMaxSizeForDirectRun = 900 # File size in MB. Any WGS file smaller than this
							 # does not need to made into smaller WGS files.

c_iReadsForFile = 5000000 # Number of WGS reads to process at a time

################################################################################
# Args

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n \
This program takes a set of protein family markers and wgs file as input, \
and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='strMarkers',
help='Enter the path and name of the genes of interest file (protein seqs).')
parser.add_argument('--wgs', type=str, dest='strWGS',
help='Enter the path and name of the WGS file (nucleotide reads).')

#Output
parser.add_argument('--results', type=str, dest='strResults', default = "results.txt",
help='Enter the name of the results file.')
parser.add_argument('--SBhits', type=str, dest='strHits',
help='ShortBRED will print the hits it considers positives to this file.', default="")
parser.add_argument('--blastout', type=str, dest='strBlast', default="",
help='Enter the name of the blast-formatted output file from USEARCH.')
parser.add_argument('--marker_results', type=str, dest='strMarkerResults', default="markers.tab",
help='Enter the name of the output for marker level results.')


parser.add_argument('--tmp', type=str, dest='strTmp', default ="",help='Enter the path and name of the tmp directory.')

#Parameters - Matching Settings
parser.add_argument('--id', type=float, dest='dID', help='Enter the percent identity for the match', default = .95)
parser.add_argument('--pctlength', type=float, dest='dAlnLength', help='Enter the minimum alignment length. The default is 20', default = 0.95)
parser.add_argument('--minreadBP', type=float, dest='iMinReadBP', help='Enter the lower bound for read lengths that shortbred witll process', default = 90)
parser.add_argument('--avgreadBP', type=float, dest='iAvgReadBP', help='Enter the average read length.', default = 100)
#parser.add_argument('--tmid', type=float, dest='dTMID', help='Enter the percent identity for a TM match', default = .95)
#parser.add_argument('--qmid', type=float, dest='dQMID', help='Enter the percent identity for a QM match', default = .95)
#parser.add_argument('--alnTM', type=int, dest='iAlnMax', help='Enter a bound for TM alignments, such that aln must be>= min(markerlength,alnTM)', default = 20)

#Parameters - Matching Various
parser.add_argument('--bz2', type=bool, dest='fbz2file', help='Set to True if using a tar.bz2 file', default = False)
parser.add_argument('--threads', type=int, dest='iThreads', help='Enter the number of CPUs available for USEARCH.', default=1)
parser.add_argument('--notmarkers', type=str, dest='strCentroids',default="N", help='This flag is used when testing centroids for evaluation purposes.')
parser.add_argument('--small', type=bool, dest='bSmall',default=False, help='This flag is used to indicate the input file is small enough for USEARCH.')

#parser.add_argument('--length', type=int, dest='iLength', help='Enter the minimum length of the markers.')

args = parser.parse_args()

if (args.strMarkers == "" or args.strWGS==""):
	parser.print_help( )
	raise Exception( "Command line arguments incorrect, must provide:\n" +
		"\t--markers AND --wgs, \n")


################################################################################
#Make temp directory
dirTmp = args.strTmp
if(dirTmp==""):
	# dirTmp gets a pid and timestamp. (This is to avoid overwriting files if
	# someone launches multiple instances of the program.)
    dirTmp = ("tmp" + str(os.getpid()) + '%.0f' % round((time.time()*1000), 1))

dirTmp = src.check_create_dir( dirTmp )

# Assign file names
if args.strHits != "":
	strHitsFile = args.strHits
else:
	strHitsFile = ( dirTmp + os.sep + "SBhits.txt" )


strMarkerResults = args.strMarkerResults
if strMarkerResults == "":
	strMarkerResults = dirTmp + os.sep + "markers.tab"


##############################################################################
# Log the parameters

with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "w") as log:
	log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n SEARCH PARAMETERS \n")
	log.write("Match ID:" + str(args.dID) + "\n")
	log.write("Pct Length for Match:" + str(args.dAlnLength) + "\n")
	if args.strCentroids=="Y":
		log.write("Sequences: Centroids\n")
	else:
		log.write("Sequences: Markers\n")

##############################################################################
#Initialize Dictionaries

dictBLAST = {}
dictMarkerLen = {}
dictMarkerLenAll = {}
dictMarkerCount = {}
dictHitsForMarker = {}

###############################################################################
#Sum up the marker lengths by family, put them in a dictionary.

for seq in SeqIO.parse(args.strMarkers, "fasta"):
	#For ShortBRED Markers...
	if args.strCentroids=="Y":
		strStub = seq.id
	#For Centroids...
	else:
		mtchStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',seq.id)
		strStub = mtchStub.group(1)

	dictMarkerLenAll[strStub] = len(seq) + dictMarkerLenAll.get(strStub,0)
	dictMarkerCount[strStub] = dictMarkerCount.get(strStub,0) + 1
	dictHitsForMarker[seq.id] = 0
 	dictMarkerLen[seq.id] = len(seq)

###############################################################################
# Make the USEARCH Database, run USEARCH, Print the ShortBRED Counts

#Make a database from the markers
strDBName = str(dirTmp) + os.sep + os.path.basename(str(args.strMarkers)) + ".udb"
sq.MakedbUSEARCH (args.strMarkers, strDBName)

#FIX THIS ONE!
if (args.strBlast == ""):
	strBlast = str(dirTmp) + os.sep + "full_results.tab"
else:
	strBlast = args.strBlast

iTotalReadCount = 0
dAvgReadLength  = 0.0

dFileInMB = round(os.path.getsize(args.strWGS)/1048576.0,1)

sys.stderr.write("The WGS file size is " + str(dFileInMB) + " MB \n")

#Check the extension on the WGS fasta file, choose appropriate extraction method.
if args.strWGS.find(".tar.bz2") > -1:
	strExtractMethod = 'r:bz2'
elif args.strWGS.find(".tar.gz") > -1:
	strExtractMethod = 'r:gz'
elif args.strWGS.find(".gz") > -1:
	strExtractMethod = 'gz'
elif args.strWGS.find(".bz2") > -1:
	strExtractMethod = 'bz2'
else:
	strExtractMethod = ""

# If the file is small and does not need to be extracted, just process directly in USEARCH.
if (dFileInMB < c_iMaxSizeForDirectRun and strExtractMethod== ""):
	args.bSmall = True
	sq.RunUSEARCH(strMarkers=args.strMarkers, strWGS=args.strWGS,strDB=strDBName, strBlastOut = strBlast,iThreads=args.iThreads,dID=args.dID, dirTmp=dirTmp )
	iMin = 999 #Can be any integer. Just a dummy to initialize iMin before calculations begin.
	for seq in SeqIO.parse(args.strWGS, "fasta"):
		iTotalReadCount+=1
		dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)
		iMin = min(iMin,len(seq))
	sq.StoreHitCounts(strBlastOut = strBlast,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
	dictHitCounts=dictBLAST,dID=args.dID,strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength,iMinReadAA=int(math.floor(args.iMinReadBP/3)),
	iAvgReadAA=int(math.floor(args.iAvgReadBP/3)))



# Otherwise, break up the large file into several small fasta files, process each one.
else:
	iReadsInSmallFile = 0
	iFileCount = 1

	strFASTAName = str(dirTmp) + os.sep + 'fasta.fna'




	#Get the list of files inside the tar file
	if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
		tarWGS = tarfile.open(args.strWGS,strExtractMethod)
		astrFileList = tarWGS.getnames()
	else:
		astrFileList = []
		astrFileList.append(args.strWGS)

	# Open each file in astrFileList with the appropriate method, pass to stream WGS
	for strFile in astrFileList:
			if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
				sys.stderr.write("Unpacking tar file... This often takes several minutes. ")
				streamWGS = tarWGS.extractfile(strFile)
			elif strExtractMethod== 'gz':
				sys.stderr.write("Unpacking gz file... This may several minutes. ")
				streamWGS = gzip.open(strFile, 'rb')
			elif strExtractMethod== 'bz2':
				sys.stderr.write("Unpacking bz2 file... This may several minutes. ")
				streamWGS =  bz2.BZ2File(strFile)
			else:
				streamWGS = open(strFile,'r')

			# While there are seqs in the current file (could be one of many if
			# the original input is tar,bz,etc.) write seqs to a small, temporary,
			# fasta file. Once the file has c_iReadsForFile, process it.


			if strFile.find("fastq") > -1:
				strFormat = "fastq"
			elif strFile.find("fasta") > -1:
				strFormat = "fasta"
			else:
				strFormat = "unknown"

			if strFormat != "unknown":
				fileFASTA = open(strFASTAName, 'w')
				iMin = 999
				for seq in SeqIO.parse(streamWGS, strFormat):
					SeqIO.write(seq,fileFASTA,"fasta")
					iReadsInSmallFile+=1
					iTotalReadCount+=1

	    			# Have a running average of the read length. This covers all of the reads in the original input file.
					dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)

					# Track minimum...
					iMin = min(len(seq),iMin)


					if (iReadsInSmallFile>=c_iReadsForFile):
						fileFASTA.close()

						#Run Usearch, store results
						strOutputName = str(dirTmp) + os.sep + "wgsout_" + str(iFileCount).zfill(2) + ".out"
						sq.RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName,dirTmp=dirTmp,iThreads=args.iThreads,dID=args.dID )
						sq.StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile,dictHitsForMarker=dictHitsForMarker, dictMarkerLen=dictMarkerLen,
						dictHitCounts=dictBLAST,dID=args.dID,strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength,iMinReadAA=int(math.floor(args.iMinReadBP/3)),
						iAvgReadAA=int(math.floor(args.iAvgReadBP/3)))

						#Reset count, make new file
						iReadsInSmallFile = 0
						iFileCount+=1
						fileFASTA = open(strFASTAName, 'w')


				if(iReadsInSmallFile>0):
					fileFASTA.close()
					#Run Usearch, store results
					strOutputName = str(dirTmp) + os.sep + "wgsout_" + str(iFileCount).zfill(2) + ".out"
					sq.RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName,dirTmp=dirTmp,iThreads=args.iThreads,dID=args.dID )
					sq.StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
					dictHitCounts=dictBLAST,dID=args.dID,strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength,iMinReadAA=int(math.floor(args.iMinReadBP/3)),
					iAvgReadAA=int(math.floor(args.iAvgReadBP/3)))

sq.CalculateCounts(strResults = args.strResults, strMarkerResults=strMarkerResults,dictHitCounts=dictBLAST,
dictMarkerLenAll=dictMarkerLenAll,dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
dReadLength = float(args.iAvgReadBP), iWGSReads = iTotalReadCount, strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength)

# Add final details to log
with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "a") as log:
	log.write("Total Reads Processed: " + str(iTotalReadCount) + "\n")
	log.write("Average Read Length: " + str(dAvgReadLength) + "\n")
	log.write("Min Read Length: " + str(iMin) + "\n")

if args.bSmall == False:
	#Delete the small, temp fasta file.
	os.remove(strFASTAName)

sys.stderr.write("Processing complete. \n")