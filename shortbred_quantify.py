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

c_iReadsForFile = 7000000 # Number of WGS reads to process at a time

################################################################################
# Args

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n \
This program takes a set of protein family markers and wgs file as input, \
and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='strMarkers',
help='Enter the path and name of the genes of interest file (protein seqs).')
parser.add_argument('--wgs', type=str, dest='strWGS',nargs='+',
help='Enter the path and name of the WGS file (nucleotide reads).')


#Output
parser.add_argument('--results', type=str, dest='strResults', default = "results.txt",
help='Enter the name of the results file.')
parser.add_argument('--SBhits', type=str, dest='strHits',
help='ShortBRED will print the hits it considers positives to this file.', default="")
parser.add_argument('--blastout', type=str, dest='strBlast', default="",
help='Enter the name of the blast-formatted output file from USEARCH.')
parser.add_argument('--marker_results', type=str, dest='strMarkerResults', default="",
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

strLog = str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log")
with open(strLog, "w") as log:
	log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n SEARCH PARAMETERS \n")
	log.write("Match ID:" + str(args.dID) + "\n")
	log.write("Pct Length for Match:" + str(args.dAlnLength) + "\n")
	if args.strCentroids=="Y":
		log.write("Sequences: Centroids\n")
	else:
		log.write("Sequences: Markers\n")

##############################################################################
#Initialize Dictionaries, Some Output Files

dictBLAST = {}
dictMarkerLen = {}
dictMarkerLenAll = {}
dictMarkerCount = {}
dictHitsForMarker = {}

#FIX THIS ONE!
if (args.strBlast == ""):
	strBlast = str(dirTmp) + os.sep + "full_results.tab"
else:
	strBlast = args.strBlast

###############################################################################
#Step 1: Prepare markers.
# Sum up the marker lengths by family, put them in a dictionary.
# Make them into a USEARCH database.

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

#Make a database from the markers
strDBName = str(dirTmp) + os.sep + os.path.basename(str(args.strMarkers)) + ".udb"
sq.MakedbUSEARCH (args.strMarkers, strDBName)

##################################################################################
#Step 2: Get information on WGS file(s), put it into aaFileInfo.

"""
aaFileInfo is array of string arrays, each with details on the file so ShortBRED
knows how to process it efficiently. Each line has the format:
	[filename, format, "large" or "small", extract method, and corresponding tarfile (if needed)]

An example:
    ['SRS011397/SRS011397.denovo_duplicates_marked.trimmed.1.fastq', 'fastq', 'large', 'r:bz2', '/n/CHB/data/hmp/wgs/samplesfqs/SRS011397.tar.bz2']
"""
astrWGS = args.strWGS

sys.stderr.write( "\nList of files in WGS set:")
for strWGS in astrWGS:
	sys.stderr.write( strWGS + "\n")

aaWGSInfo = []

for strWGS in astrWGS:
	strExtractMethod= sq.CheckExtract(strWGS)

	# If tar file, get details on members, and note corresponding tarfile
	# Remember that a tarfile has a header block, and then data blocks
	if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
		tarWGS = tarfile.open(strWGS,strExtractMethod)
 		atarinfoFiles = tarWGS.getmembers() #getmembers() returns tarInfo objects
		tarWGS.close()

		for tarinfoFile in atarinfoFiles:
			if tarinfoFile.isfile(): # This condition confirms that it is a file, not a header.
				strFormat = sq.CheckFormat(tarinfoFile.name)
				strSize = sq.CheckSize(tarinfoFile.size, c_iMaxSizeForDirectRun)
				astrFileInfo = [tarinfoFile.name, strFormat, strSize,strExtractMethod, strWGS ]
				aaWGSInfo.append(astrFileInfo)
    # Otherwise, get file details directly
	else:
		strFormat = sq.CheckFormat(strWGS)
		dFileInMB = round(os.path.getsize(strWGS)/1048576.0,1)
		if dFileInMB < c_iMaxSizeForDirectRun:
			strSize = "small"
		else:
			strSize = "large"
		astrFileInfo = [strWGS, strFormat, strSize,strExtractMethod, "no_tar" ]
		aaWGSInfo.append(astrFileInfo)


sys.stderr.write( "\nList of files in WGS set (after unpacking tarfiles):")
for astrWGS in aaWGSInfo:
	sys.stderr.write( astrWGS[0]+" ")

sys.stderr.write("\n\n")
##################################################################################
# Step 3: Call USEARCH on each WGS file, (break into smaller files if needed), store hit counts.

# Initialize values for the sample
iTotalReadCount = 0
dAvgReadLength  = 0.0
iMin = 999 #Can be any integer. Just a dummy to initialize iMin before calculations begin.
iWGSFileCount = 1

with open(strLog, "a") as log:
	log.write('\t'.join(["# FileName","size","format","extract method","tar file (if part of one)"]) + '\n')
	log.write("Reads processed" + "\n")

for astrFileInfo in aaWGSInfo:
	strWGS,strFormat,strSize,strExtractMethod,strMainTar = astrFileInfo
	with open(strLog, "a") as log:
		log.write(str(iWGSFileCount) + ": " + '\t'.join(astrFileInfo) + '\n')
	iWGSReads = 0
	sys.stderr.write( "Working on file " + str(iWGSFileCount) + " of " + str(len(aaWGSInfo)))

	#If it's a small fasta file, just give it to USEARCH directly.
	if strSize=="small" and strFormat=="fasta":
		sq.RunUSEARCH(strMarkers=args.strMarkers, strWGS=strWGS,strDB=strDBName, strBlastOut = strBlast,iThreads=args.iThreads,dID=args.dID, dirTmp=dirTmp )
		sq.StoreHitCounts(strBlastOut = strBlast,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
			dictHitCounts=dictBLAST,dID=args.dID,strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength,iMinReadAA=int(math.floor(args.iMinReadBP/3)),
			iAvgReadAA=int(math.floor(args.iAvgReadBP/3)))


		for seq in SeqIO.parse(strWGS, "fasta"):
			iWGSReads+=1
			iTotalReadCount+=1
			dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)
			iMin = min(iMin,len(seq))


		"""
		#Skip the file if the format is unknown.
		elif strFormat == "unknown":
			sys.stderr.write("WARNING: Skipped file with unknown format: " + strWGS + "\n")
			with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "a") as log:
				log.write("WARNING: Skipped file with unknown format: " + strWGS + "\n")
		"""

	#Otherwise, convert the file as needed to into small fasta files. Call USEARCH and store the counts for each small file.
	else:
		iReadsInSmallFile = 0
		iFileCount = 1

		strFASTAName = str(dirTmp) + os.sep + 'fasta.fna'

		#Unpack file with appropriate extract method
		if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
			sys.stderr.write("Unpacking tar file... this often takes several minutes. ")
			tarWGS = tarfile.open(strMainTar,strExtractMethod)
			streamWGS = tarWGS.extractfile(strWGS)
		elif strExtractMethod== 'gz':
			sys.stderr.write("Unpacking gz file... this may take several minutes. ")
			streamWGS = gzip.open(strWGS, 'rb')
		elif strExtractMethod== 'bz2':
			sys.stderr.write("Unpacking bz2 file... this may take several minutes. ")
			streamWGS =  bz2.BZ2File(strWGS)
		else:
			streamWGS = open(strWGS,'r')

		#Open file for writing
		fileFASTA = open(strFASTAName, 'w')

		#TEMPORARY TEST!
		if strFormat=="unknown":
			strFormat="fastq"
		if streamWGS==None:
			with open(strLog, "a") as log:
				log.write("File was empty." + '\n')

		#Start the main loop to get everything in streamWGS -> small fasta file -> counted and stored
		for seq in SeqIO.parse(streamWGS, strFormat):
			SeqIO.write(seq,fileFASTA,"fasta")
			iReadsInSmallFile+=1
			iTotalReadCount+=1
			iWGSReads+=1

			# Have a running average of the read length. This covers all of the reads in the original input file.
			dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)
			iMin = min(len(seq),iMin)

			#Close the temp fasta file once it has enough reads.
			if (iReadsInSmallFile>=c_iReadsForFile):
				fileFASTA.close()

				#Run Usearch, store results
				strOutputName = str(dirTmp) + os.sep + "wgs_" + str(iWGSFileCount).zfill(2) + "out_" + str(iFileCount).zfill(2) + ".out"
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
			strOutputName = str(dirTmp) + os.sep + "wgs_" + str(iWGSFileCount).zfill(2) + "out_" + str(iFileCount).zfill(2) + ".out"
			sq.RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName,dirTmp=dirTmp,iThreads=args.iThreads,dID=args.dID )
			sq.StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
			dictHitCounts=dictBLAST,dID=args.dID,strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength,iMinReadAA=int(math.floor(args.iMinReadBP/3)),
			iAvgReadAA=int(math.floor(args.iAvgReadBP/3)))

	with open(strLog, "a") as log:
		log.write(str(iWGSReads) + '\n')

	iWGSFileCount += 1
	if (strFormat != "fasta" or strSize != "small"):
		streamWGS.close()
	#Close the tarfile if you had one open.
	if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
		tarWGS.close()
##################################################################################
# Step 4: Calculate ShortBRED Counts, print results, print log info.

sq.CalculateCounts(strResults = args.strResults, strMarkerResults=strMarkerResults,dictHitCounts=dictBLAST,
dictMarkerLenAll=dictMarkerLenAll,dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,
dReadLength = float(args.iAvgReadBP), iWGSReads = iTotalReadCount, strCentCheck=args.strCentroids,dAlnLength=args.dAlnLength)

# Add final details to log
with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "a") as log:
	log.write("Total Reads Processed: " + str(iTotalReadCount) + "\n")
	log.write("Average Read Length: " + str(dAvgReadLength) + "\n")
	log.write("Min Read Length: " + str(iMin) + "\n")

if strSize != "small":
	#Delete the small, temp fasta file.
	os.remove(strFASTAName)

sys.stderr.write("Processing complete. \n")
