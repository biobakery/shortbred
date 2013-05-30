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
import time

import src
import numpy

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

c_strUSEARCH	= "usearch"

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers and wgs file as input, \
and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='strMarkers', help='Enter the path and name of the genes of interest file (protein seqs).')
parser.add_argument('--wgs', type=str, dest='strWGS', help='Enter the path and name of the WGS file (nucleotide reads).')

#Output
parser.add_argument('--results', type=str, dest='strResults', default = "results.txt",help='Enter the name of the results file.')
parser.add_argument('--SBhits', type=str, dest='strHits', help='ShortBRED will print the hits it considers positives to this file.', default="")
parser.add_argument('--blastout', type=str, dest='strBlast', default="out.blast",help='Enter the name of the blast-formatted output file from USEARCH.')
parser.add_argument('--marker_results', type=str, dest='strMarkerResults', default="markers.tab",help='Enter the name of the output for marker level results.')


parser.add_argument('--tmp', type=str, dest='strTmp', default ="",help='Enter the path and name of the tmp directory.')

#Parameters - Matching Settings
parser.add_argument('--id', type=float, dest='dID', help='Enter the percent identity for the match', default = .95)
parser.add_argument('--tmid', type=float, dest='dTMID', help='Enter the percent identity for a TM match', default = .95)
parser.add_argument('--qmid', type=float, dest='dQMID', help='Enter the percent identity for a QM match', default = .95)
parser.add_argument('--alnlength', type=int, dest='iAlnLength', help='Enter the minimum alignment length. The default is 20', default = 20)
parser.add_argument('--alnTM', type=int, dest='iAlnMax', help='Enter a bound for TM alignments, such that aln must be>= min(markerlength,alnTM)', default = 20)

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

if args.strHits != "":
	strHitsFile = args.strHits
else:
	strHitsFile = ( dirTmp + os.sep + "SBhits.txt" )

###############################################################################




def RunUSEARCH ( strMarkers, strWGS,strBlastOut, strDB):

	subprocess.check_call(["time","-o", strMarkers + ".time",
		c_strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(args.dID),"--blast6out", strBlastOut,
		"--threads", str(args.iThreads)])

def StoreHitCounts(strBlastOut,strValidHits,dictHitsForMarker,dictMarkerLen,dictHitCounts):
	#strBlastOut - BLAST-formatted output from USEARCH
	#strValidHits - File of BLAST hits that meet ShortBRED's ID and Length criteria. Mainly used for evaluation/debugging.
	#dictMarkerLen - Contains each marker/centroid length
	#dictHitCounts - Contains each family's hit count


	csvwHits = csv.writer( open(strValidHits,'a'), csv.excel_tab )

	#Go through the usearch output, for each prot family, record the number of valid hits
	for aLine in csv.reader( open(strBlastOut), csv.excel_tab ):

		#Pick appropriate ID and Length criteria, based on whether marker is a TM or QM
		mtchTM = re.search(r'_TM',aLine[1])
		#Changed so that QM's and TM's are treated the same.
		"""
		if (mtchTM):
			dID = args.dTMID
			iAln = min(dictMarkerLen[aLine[1]] ,args.iAlnMax)
		else:
			dID = args.dQMID
			iAln = args.iAlnLength
		"""
		dID = args.dTMID
		iAln = min(dictMarkerLen[aLine[1]] ,args.iAlnMax)

		#If using centroids (Typically only used for evaluation purposes.)....
		if args.strCentroids=="Y":
			strProtFamily = aLine[1]

			if (int(aLine[3])>= iAln):
					dictHitCounts[strProtFamily] = dictHitCounts.setdefault(strProtFamily,0) + 1
					dictHitsForMarker[strProtFamily] = dictHitsForMarker.setdefault(strProtFamily,0) + 1
					csvwHits.writerow( aLine )

		#If using ShortBRED Markers (and not centroids)...
		else:
			#Get the Family Name
			mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])
			strProtFamily = mtchProtStub.group(1)

			#If hit satisfies criteria, add it to dictHitCounts's count of hits for that family, write the result to fileHits.
			if (int(aLine[3])>= iAln and (float(aLine[2])/100) >= dID):
				#Add 1 to count of hits for that marker
				dictHitsForMarker[aLine[1]] = dictHitsForMarker.setdefault(aLine[1],0) + 1

                #Add 1 to count of hits for that family
				dictHitCounts[strProtFamily] = dictHitCounts.setdefault(strProtFamily,0) +1

				csvwHits.writerow( aLine )


def ProcessHitData(atupHits):
	strMarkerFile = args.strMarkerResults
	if strMarkerFile == "":
		strMarkerFile = dirTmp + os.sep + "markers.tab"
	with open(strMarkerFile, 'w') as csvfileMarker:
		csvwMarkerResults = csv.writer( csvfileMarker, csv.excel_tab )
		csvwMarkerResults.writerow(["Family","Marker","Normalized Count","Hits","MarkerLength","ReadLength"])

	strFamFile = args.strResults
	if strFamFile == "":
		strFamFile = dirTmp + os.sep + "families.tab"
	with open(strFamFile, 'w') as csvfileFam:
		csvwFamResults = csv.writer( csvfileFam, csv.excel_tab )
		csvwFamResults.writerow(["Family","Count","Hits","TotMarkerLength"])

	# Sort them by Family Name
	atupHits.sort(key=lambda x: x[0])



	strCurFam = ""
	atupCurFamData = []


	for tupRow in atupHits:
		strFam = tupRow[0]
		if strFam != strCurFam:
			# Print results, start a new array for this family.
			if strCurFam!="":
				PrintStats(atupCurFamData, strMarkerFile,strFamFile)
			strCurFam = strFam
			atupCurFamData = []
			atupCurFamData.append(tupRow)
		else:
			# Add to the current array.
			atupCurFamData.append(tupRow)


	PrintStats(atupCurFamData, strMarkerFile,strFamFile)
	return

def PrintStats(atupCurFamData, strMarkerFile, strFamFile):


	# We want two files:
	# 1) stats by family
	# 2) stats for each marker (sorted by family and marker.)
	#for now...
	#Sort by the marker
	atupCurFamData.sort(key=lambda x: x[1])


	with open(strMarkerFile, 'a') as csvfileMarker:
		csvwMarkerResults = csv.writer( csvfileMarker, csv.excel_tab )
		for tupRow in atupCurFamData:
			csvwMarkerResults.writerow(tupRow)
			strName = tupRow[0]

	#Zip the tuples so that we perform operations on the columns.
	atupZipped = zip(*atupCurFamData)


	# Family Stats
	try:
		dMedian = numpy.median(list(atupZipped[2]))
		iHits = sum(list(atupZipped[3]))
		iMarkerLength = sum(list(atupZipped[4]))
	except:
         sys.stderr.write("Problem with results for set:",str(atupZipped))


	# Print out the family results
	with open(strFamFile, 'a') as csvfile:
		csvwFamResults = csv.writer( csvfile, csv.excel_tab )
		csvwFamResults.writerow([strName,dMedian,iHits,iMarkerLength])

	return

def PrintResults(strResults,dictHitCounts, dictHitsForMarker, dictMarkerLenAll,dictMarkerLen,dReadLength,iWGSReads):
	#strResults - Name of text file with final ShortBRED Counts
	#strBlastOut - BLAST-formatted output from USEARCH
	#strValidHits - File of BLAST hits that meet ShortBRED's ID and Length criteria. Mainly used for evaluation/debugging.
	#dictMarkerLenAll - Contains the sum of marker lengths for all markers in a family
	#dictMarkerLen - Contains each marker/centroid length

	atupMarkerCounts = []

	#Print Name, Normalized Count, Hit Count, Marker Length to std out
	#csvwResults = csv.writer( open(strResults,'w'), csv.excel_tab )
	#csvwResults.writerow(["Marker","Normalized Count","Hits","MarkerLength","ReadLength"])

	for strMarker in dictHitsForMarker.keys():
		# Switching to Method 1
		"""
  			If marker length < average read length:
                  Count = Hits /  [AvgReadLength - MarkerLenInNucs) / AvgReadLength]
             else:
                  Count = Hits /  [MarkerLenInNucs - AvgReadLength) / AvgReadLength]
		"""
		iMarkerNucs = dictMarkerLen[strMarker]*3
		iAlnLength = args.iAlnMax*3
		iHits = dictHitsForMarker[strMarker]

		# Consider the problem as fitting shorter sequence into the longer sequence.
		# We add 1 for the special case when the marker is as long as the read.
		dCount = iHits/ ( (abs(dReadLength - iMarkerNucs)+1) / float(dReadLength))

		# Normalize for metagenome depth
		#dCount = (dCount / (iWGSReads))*1000
		dCount = dCount * 1000 / (iWGSReads / 1e9)

		if args.strCentroids=="Y":
			strProtFamily = strMarker
		else:
			mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
			strProtFamily = mtchProtStub.group(1)


		tupCount = (strProtFamily,strMarker, dCount,dictHitsForMarker[strMarker],dictMarkerLen[strMarker],dReadLength)
		atupMarkerCounts.append(tupCount)

        ProcessHitData(atupMarkerCounts)



	"""
	OLD METHOD
	for strProt in dictHitCounts.keys():
		dShortBREDCount = (float(dictHitCounts[strProt])/(dictMarkerLenAll[strProt]) - (dictMarkerCount * (dAvgReadLength/3)) / (float(iWGSReads) /dReadLength))
		csvwResults.writerow( [strProt, dShortBREDCount,
			dictHitCounts[strProt], dictMarkerLenAll[strProt]] )
	"""
"""
reads hitting markers / [total length of markers - (avg read length * # of markers longer than avg read length)]
"""

##############################################################################
# Log the parameters

with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "w") as log:
	log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n SEARCH PARAMETERS \n")
	log.write("Match ID:" + str(args.dID) + "\n")
	log.write("Alignment Length:" + str(args.iAlnLength) + "\n")
	log.write("TM id:" + str(args.dTMID) + "\n")
	log.write("QM id:" + str(args.dQMID) + "\n")
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

p = subprocess.check_call([c_strUSEARCH, "--makeudb_usearch", args.strMarkers,
	"--output", strDBName])

strBlast = args.strBlast

if (args.bSmall == True):
	iTotalReadCount = 0
	dAvgReadLength = 0.0
	RunUSEARCH(strMarkers=args.strMarkers, strWGS=args.strWGS,strDB=strDBName, strBlastOut = strBlast )
	StoreHitCounts(strBlastOut = strBlast,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)

	for seq in SeqIO.parse(args.strWGS, "fasta"):
		iTotalReadCount+=1
		dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)


else:
	#Check the extension on the WGS fasta file
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

	iReadsForFile = 5000000
	iCount = 0
	iTotalReadCount = 0
	dAvgReadLength = 0.0

	iFileCount = 1
	strFASTAName = str(dirTmp) + os.sep + 'fasta.fna'
	fileFASTA = open(strFASTAName, 'w')

	astrFileList = []

	#Get the list of files inside the tar file
	if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
		tarWGS = tarfile.open(args.strWGS,strExtractMethod)
		astrFileList = tarWGS.getnames()
		#print(tarWGS.getnames)
	else:
		astrFileList.append(args.strWGS)

	#print "Compression: ", strExtractMethod
	#print astrFileList


	#Note: I would like to add code here along the lines of
	#if(strExtractMethod=="" and OneFile and FileISSmall):
	#   Pass to usearch without processing.

	#Open each one with the appropriate method



	for strFile in astrFileList:
			if (strExtractMethod== 'r:bz2' or strExtractMethod=='r:gz'):
				sys.stderr.write("Unzipping tar file... This often takes several minutes. ")
				streamWGS = tarWGS.extractfile(strFile)
			elif strExtractMethod== 'gz':
				sys.stderr.write("Unzipping gz file... This may several minutes. ")
				streamWGS = gzip.open(strFile, 'rb')
			elif strExtractMethod== 'bz2':
				sys.stderr.write("Unzipping bz2 file... This may several minutes. ")
				streamWGS =  bz2.BZ2File(strFile)
			else:
				streamWGS = open(strFile,'r')


			for seq in SeqIO.parse(streamWGS, "fasta"):
				SeqIO.write(seq,fileFASTA,"fasta")
				iCount+=1
				iTotalReadCount+=1

				# This tracks a running average of the read length. Illumina is typically constant,
				# but 454 reads vary in size.

				dAvgReadLength = ((dAvgReadLength * (iTotalReadCount-1)) + len(seq))/float(iTotalReadCount)



				#print iCount
				if (iCount>=iReadsForFile):
					fileFASTA.close()

					#Run Usearch, store results
					strOutputName = str(dirTmp) + os.sep + "wgsout_" + str(iFileCount).zfill(2) + ".out"
					RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName )
					StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile,dictHitsForMarker=dictHitsForMarker, dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)

					#Reset count, make new file
					iCount = 0
					iFileCount+=1
					fileFASTA = open(strFASTAName, 'w')


			if(iCount>0):
				fileFASTA.close()
				#Run Usearch, store results
				strOutputName = str(dirTmp) + os.sep + "wgsout_" + str(iFileCount).zfill(2) + ".out"
				RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName )
				StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile, dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)

PrintResults(strResults = args.strResults, dictHitCounts=dictBLAST, dictMarkerLenAll=dictMarkerLenAll,dictHitsForMarker=dictHitsForMarker,dictMarkerLen=dictMarkerLen, dReadLength = dAvgReadLength, iWGSReads = iTotalReadCount)
with open(str(dirTmp + os.sep + os.path.basename(args.strMarkers)+ ".log"), "a") as log:
	log.write("Total Reads Processed: " + str(iTotalReadCount) + "\n")
	log.write("Average Read Length: " + str(dAvgReadLength) + "\n")

if args.bSmall == False:
	os.remove(strFASTAName)
"""
***** copy to tmp fasta
***** if count exceeds iReadsForFile
****** stop
***** close tmp fasta
***** run usearch
***** store counts
**** print (aggregate) results
"""

#print iCount

"""
    	for seq in SeqIO.parse(tarWGS.extractfile(wgsFile), "fasta"):
			if (iCount< iReadsForFile):
				SeqIO.write(seq, fileFASTA, "fasta")
				iCount+=1
			else:
				fileFASTA.close()
				print iCount
				print "Making a new fasta file..."
				strOutputName = str(dirTmp) + os.sep + "wgsout_" + str(iFileCount).zfill(2) + ".out"
				RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName )
				StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile, dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)

				fileFASTA = open(strFASTAName, 'w')
				iFileCount +=1
				SeqIO.write(seq, fileFASTA, "fasta")
				iCount = 1
	fileFASTA.close()
	RunUSEARCH(strMarkers=args.strMarkers, strWGS=strFASTAName,strDB=strDBName, strBlastOut = strOutputName )
	StoreHitCounts(strBlastOut = strOutputName,strValidHits=strHitsFile, dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)
	PrintResults(strResults = args.strResults, dictHitCounts=dictBLAST, dictMarkerLenAll=dictMarkerLenAll,dictMarkerLen=dictMarkerLen)
	tarWGS.close()




#If WGS is in a fasta file...
if args.fbz2file == False:
	RunUSEARCH(strMarkers=args.strMarkers, strWGS=args.strWGS,strDB=strDBName, strBlastOut = args.strBlast )
	StoreHitCounts(strBlastOut = args.strBlast,strValidHits=strHitsFile, dictMarkerLen=dictMarkerLen,dictHitCounts=dictBLAST)
	PrintResults(strResults = args.strResults, dictHitCounts=dictBLAST, dictMarkerLenAll=dictMarkerLenAll,dictMarkerLen=dictMarkerLen)

#If WGS is in a tar file...
else:




	tarWGS = tarfile.open(args.strWGS,'r:bz2')
	for wgsFile in tarWGS.getnames():
"""


