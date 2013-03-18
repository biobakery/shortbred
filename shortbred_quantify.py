#!/usr/bin/env python

import sys
import argparse
import subprocess
import csv
import re
import os
import datetime
import shutil
import tarfile

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers and wgs file as input, and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='strMarkers', help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--wgs', type=str, dest='strWGS', help='Enter the path and name of the genes of interest file (proteins).')

#Output
parser.add_argument('--results', type=str, dest='strResults', default = "results.txt",help='Enter the name of the results file.')
parser.add_argument('--SBhits', type=str, dest='strHits', help='ShortBRED will print the hits it considers positives to this file.', default="")
parser.add_argument('--blastout', type=str, dest='strBlast', default="out.blast",help='Enter the path and name of the blastoutput.')

#Parameters
parser.add_argument('--bz2', type=bool, dest='fbz2file', help='Set to true if using a bz2 file', default = False)
parser.add_argument('--id', type=float, dest='dID', help='Enter the percent identity for the match', default = .95)
parser.add_argument('--tmid', type=float, dest='dTMID', help='Enter the percent identity for a TM match', default = .95)
parser.add_argument('--qmid', type=float, dest='dQMID', help='Enter the percent identity for a QM match', default = .90)
parser.add_argument('--alnlength', type=int, dest='iAlnLength', help='Enter the minimum alignment length. The default is 20', default = 20)
parser.add_argument('--alnTM', type=int, dest='iAlnMax', help='Enter a bound for TM alignments, such that aln must be>= min(markerlength,alnTM)', default = 25)
parser.add_argument('--tmp', type=str, dest='strTmp', default =os.getcwd() +os.sep + "tmp",help='Enter the path and name of the tmp directory.')
parser.add_argument('--length', type=int, dest='iLength', help='Enter the minimum length of the markers.')
parser.add_argument('--threads', type=int, dest='iThreads', help='Enter the number of CPUs available for usearch.', default=1)
parser.add_argument('--notmarkers', type=str, dest='strNM',default="N", help='.')

#DB Note - Maybe ask Nicola how to remove usearch6 outpu
args = parser.parse_args()

################################################################################
#Make temp directory
dirTmp = args.strTmp
if not os.path.exists(dirTmp):
    os.makedirs(dirTmp)

if(args.strHits==""):
    strHitsFile = args.strTmp + os.sep + "SBhits.txt"


###############################################################################
def RunUSEARCH ( strMarkers, strWGS,strBlastOut, strDB):

	subprocess.check_call(["time","-o", strMarkers + ".time","usearch6", "--usearch_local", strWGS, "--db", strDB, "--id", str(args.dID),"--blast6out", args.strBlast,"--threads", str(args.iThreads)])
	return

def PrintResults(strResults,strBlastOut,strValidHits, dictMarkerLenAll,dictMarkerLen):
	#strResults - Name of text file with final ShortBRED Counts
	#strBlastOut - BLAST-formatted output from USEARCH
	#strValidHits - File of BLAST hits that meet ShortBRED's ID and Length criteria. Mainly used for evaluation/debugging.
	#dictMarkerLenAll - Contains the sum of marker lengths for all markers in a family
    #dictMarkerLen - Contains each marker/centroid length

	fileHits = open(strValidHits,'w')

	#Go through the usearch output, for each prot family, record the number of valid hits
	dictBLAST = {}
	for aLine in csv.reader( open(strBlastOut), csv.excel_tab ):

		#Pick appropriate ID and Length criteria, based on whether marker is a TM or QM
	    mtchTM = re.search(r'_TM.*',aLine[1])
	    if (mtchTM):
	        dID = args.dTMID
	        iAln = min(dictMarkerLen[aLine[1]] ,args.iAlnMax)
	    else:
	        dID = args.dQMID
	        iAln = args.iAlnLength

		#If using ShortBRED Markers (and not centroids)...
	    if args.strNM=="N":
			#Get the Family Name
	        mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])
	        strProtFamily = mtchProtStub.group(1)

			#If hit satisfies criteria, add it to dictBLAST's count of hits for that family, write the result to fileHits.
	        if (int(aLine[3])>= iAln and (float(aLine[2])/100.0) >= dID):
	            iCountHits = dictBLAST.setdefault(strProtFamily,0)
	            iCountHits = iCountHits+1
	            dictBLAST[strProtFamily] = iCountHits
	            fileHits.write('\t'.join(aLine) + '\n')

		#If using centroids (mainly for evaluation)....
     	else:
	        strProtFamily = aLine[1]

	        if (int(aLine[3])>= iAln):
		            iCountHits = dictBLAST.setdefault(strProtFamily,0)
		            iCountHits = iCountHits+1
		            dictBLAST[strProtFamily] = iCountHits
		            fileHits.write('\t'.join(aLine) + '\n')

	fileHits.close()

    #Print Name, Normalized Count, Hit Count, Marker Length to std out
	fileResults = open(strResults,'w')
	for strProt in dictBLAST.keys():
	    fileResults.write(strProt + "\t" + str(float(dictBLAST[strProt])/dictMarkerLenAll[strProt]) + "\t" + str(dictBLAST[strProt]) + "\t" + str(dictMarkerLenAll[strProt]) + "\n")
	fileResults.close()

##############################################################################
# Log the parameters

log = open(str(dirTmp + os.sep + args.strMarkers+ ".log"), "w")
log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n SEARCH PARAMETERS \n")
log.write("Match ID:" + str(args.dID) + "\n")
log.write("Alignment Length:" + str(args.iAlnLength) + "\n")
log.write("TM id:" + str(args.dTMID) + "\n")
log.write("QM id:" + str(args.dQMID) + "\n")
log.close()

###############################################################################
# Sum the total marker length for each family

dictMarkerLen = {}
dictMarkerLenAll = {}

#Sum up the marker lengths by family, put them in a dictionary.
for seq in SeqIO.parse(args.strMarkers, "fasta"):
	#For ShortBRED Markers...
	if args.strNM=="N":
		mtchStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',seq.id)
		strStub = mtchStub.group(1)
    #For Centroids...
	else:
		strStub = seq.id
	dictMarkerLenAll[strStub] = len(seq) + dictMarkerLenAll.get(strStub,0)
 	dictMarkerLen[seq.id] = len(seq)

###############################################################################
# Make the USEARCH Database, run USEARCH, Print the ShortBRED Counts

#Make a database from the markers
strDBName = str(dirTmp) + os.sep + str(args.strMarkers) + ".udb"
p = subprocess.check_call(["usearch6", "--makeudb_usearch", str(args.strMarkers), "--output", strDBName])

#If WGS is in a fasta file...
if( args.fbz2file==False):
	RunUSEARCH(strMarkers=args.strMarkers, strWGS=args.strWGS,strDB=strDBName, strBlastOut = args.strBlast )
	PrintResults(strResults = args.strResults, strBlastOut= args.strBlast,strValidHits=strHitsFile, dictMarkerLenAll=dictMarkerLenAll,dictMarkerLen=dictMarkerLen)

#If WGS is in a tar file...
else:
	iLinesForFile = 500000
	iCount = 0

	fileFASTA = open(str(dirTmp) + os.sep + 'fasta.fna', 'w')

	tarWGS = tarfile.open(args.strWGS,'r:bz2')
	for wgsFile in tarWGS.getnames():
		for strLine in tarWGS.extractfile(wgsFile):
			if (iCount< iLinesForFile):
				fileFASTA.write(strLine)
				iCount+=1
			else:
				fileFASTA.close()
				print iCount
				print "Making a new fasta file..."
				fileFASTA = open(str(dirTmp) + os.sep + 'fasta.fna', 'w')
				fileFASTA.write(strLine)
				iCount = 1

	fileFASTA.close()
	tarWGS.close()


