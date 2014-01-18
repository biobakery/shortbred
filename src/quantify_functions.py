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

import subprocess
import csv
import re
import sys
import math
import os


c_iAlnCentroids = 30

def MakedbUSEARCH ( strMarkers, strDBName,strUSEARCH):
	p = subprocess.check_call([strUSEARCH, "--makeudb_usearch", strMarkers,"--output", strDBName])
	return

def CheckFormat ( strFile):
	if strFile.find("fastq") > -1:
		strFormat = "fastq"
	elif strFile.find("fasta") > -1:
		strFormat = "fasta"
	elif strFile.find(".fna") > -1:
		strFormat = "fasta"
	elif strFile.find(".faa") > -1:
		strFormat = "fasta"
	else:
		strFormat = "unknown"

	return strFormat

def CheckExtract(strWGS):
	if strWGS.find(".tar.bz2") > -1:
		strExtractMethod = 'r:bz2'
	elif strWGS.find(".tar.gz") > -1:
		strExtractMethod = 'r:gz'
	elif strWGS.find(".gz") > -1:
		strExtractMethod = 'gz'
	elif strWGS.find(".bz2") > -1:
		strExtractMethod = 'bz2'
	else:
		strExtractMethod = ""

	return strExtractMethod

def CheckSize(iSize, iMax):
	dFileInMB = round(iSize/1048576.0,1)
	if dFileInMB < iMax:
		strSize = "small"
	else:
		strSize = "large"

	return strSize

def RunUSEARCH ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):

	strFields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts"

	subprocess.check_call(["time","-o", str(dirTmp) + os.sep + os.path.basename(strMarkers) + ".time",
		strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strBlastOut,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])

def RunUSEARCHGenome ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):

	strFields = "target+query+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts"

	subprocess.check_call(["time","-o", str(dirTmp) + os.sep + os.path.basename(strMarkers) + ".time",
		strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strBlastOut,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])

def Median(adValues):
	adValues.sort()
	iLen = len(adValues)
	if iLen % 2==0:
		dMedian = float(adValues[iLen/2] + adValues[(iLen/2)-1])/2.0
	else:
		dMedian = adValues[int(math.floor(iLen/2))]
	return dMedian

def StoreHitCounts(strBlastOut,strValidHits,dictHitsForMarker,dictMarkerLen,dictHitCounts,dID,strCentCheck,dAlnLength,iMinReadAA,iAvgReadAA):
# Reads in the USEARCH output (strBlastOut), marks which hits are valid (id>=dID &
# len >= min(95% of read,dictMarkerLen[Marker]) and adds to count in dictHitsForMarker[strMarker].
# Valid hits are also copied to the file in strValidHits. strCentCheck is used to flag centroids,
# and handle their counting


	with open(strValidHits, 'a') as csvfileHits:
		csvwHits = csv.writer( csvfileHits, csv.excel_tab )


		sys.stderr.write("Processing USEARCH results... \n")
		#Go through the usearch output, for each prot family, record the number of valid

		with open(strBlastOut, 'r') as csvfileBlast:
			for aLine in csv.reader( csvfileBlast, delimiter='\t' ):


				strMarker 	= aLine[1]
				dHitID		= aLine[2]
				iAlnLen     = int(aLine[3])
				iReadLenAA  = int(aLine[12])

				# A valid match must be as long as 95% of the read or the full marker.
		        # (Note that this in AA's.)


				#If using centroids (Typically only used for evaluation purposes.)....
				if strCentCheck=="Y":
					strProtFamily = strMarker

					if ( (int(iAlnLen)>= c_iAlnCentroids) and ( float(dHitID)/100) >= dID):
							dictHitCounts[strProtFamily] = dictHitCounts.setdefault(strProtFamily,0) + 1
							dictHitsForMarker[strProtFamily] = dictHitsForMarker.setdefault(strProtFamily,0) + 1
							csvwHits.writerow( aLine )

				#If using ShortBRED Markers (and not centroids)...
				else:
					iAlnMin = min(dictMarkerLen[strMarker] ,math.floor((iAvgReadAA)*dAlnLength))
					#Get the Family Name
					mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
					strProtFamily = mtchProtStub.group(1)

					#If hit satisfies criteria, add it to counts, write out data to Hits file
					if (int(iAlnLen)>= iAlnMin and (iReadLenAA >= iMinReadAA) and (float(dHitID)/100) >= dID):

						#Add 1 to count of hits for that marker, and family
						dictHitsForMarker[aLine[1]] = dictHitsForMarker.setdefault(aLine[1],0) + 1
						dictHitCounts[strProtFamily] = dictHitCounts.setdefault(strProtFamily,0) +1

						csvwHits.writerow( aLine )
	return


"""
CalculateCounts - Calculates the ShortBRED counts for each marker.
ProcessHitData - Opens the marker and family results files for writing, calls
PrintStats for each family.
PrintStats -
"""


def ProcessHitData(atupHits,strMarkerResults,strFamFile):
# Called by CalculateCounts. This function takes the set of ShortBRED marker results,
# (atupHits) and calls PrintStats to print out their results to the Marker
# (strMarkerResults) and Family results (strFamFile).
	with open(strMarkerResults, 'w') as csvfileMarker:
		csvwMarkerResults = csv.writer( csvfileMarker, csv.excel_tab )
		csvwMarkerResults.writerow(["Family","Marker","Normalized Count","Hits","MarkerLength","ReadLength","HitSpace"])

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
				PrintStats(atupCurFamData,  strMarkerFile=strMarkerResults,strFamFile=strFamFile)
			strCurFam = strFam
			atupCurFamData = []
			atupCurFamData.append(tupRow)
		else:
			# Add to the current array.
			atupCurFamData.append(tupRow)


	PrintStats(atupCurFamData, strMarkerFile=strMarkerResults,strFamFile=strFamFile)
	return

def PrintStats(atupCurFamData, strMarkerFile, strFamFile):
# Called by ProcessHitData. This function takes the set of ShortBRED marker results
# for one family(atupCurFamData), and appends their results to the Marker
# (strMarkerResults) and Family results (strFamFile).

	atupCurFamData.sort(key=lambda x: x[1])

	# Print out the marker results
	with open(strMarkerFile, 'a') as csvfileMarker:
		csvwMarkerResults = csv.writer( csvfileMarker, csv.excel_tab )
		for tupRow in atupCurFamData:
			csvwMarkerResults.writerow(tupRow)
			strName = tupRow[0]

	#sys.stderr.write("Processing "+strName+"... \n")
	#Zip the tuples so that we perform operations on the columns.
	atupZipped = zip(*atupCurFamData)


	# Family Stats
	try:
		dMedian = Median(list(atupZipped[2]))
		iHits = sum(list(atupZipped[3]))
		iMarkerLength = sum(list(atupZipped[4]))
	except:
         sys.stderr.write("Problem with results for set: " +str(atupZipped))


	# Print out the family results
	with open(strFamFile, 'a') as csvfile:
		csvwFamResults = csv.writer( csvfile, csv.excel_tab )
		csvwFamResults.writerow([strName,dMedian,iHits,iMarkerLength])

	return

def CalculateCounts(strResults,strMarkerResults, dictHitCounts, dictHitsForMarker, dictMarkerLenAll,dictMarkerLen,dReadLength,iWGSReads,strCentCheck,
dAlnLength,strFile):
	#strResults - Name of text file with final ShortBRED Counts
	#strBlastOut - BLAST-formatted output from USEARCH
	#strValidHits - File of BLAST hits that meet ShortBRED's ID and Length criteria. Mainly used for evaluation/debugging.
	#dictMarkerLenAll - Contains the sum of marker lengths for all markers in a family
	#dictMarkerLen - Contains each marker/centroid length

	atupMarkerCounts = []

	#Print Name, Normalized Count, Hit Count, Marker Length to std out
	#csvwResults = csv.writer( open(strResults,'w'), csv.excel_tab )
	#csvwResults.writerow(["Marker","Normalized Count","Hits","MarkerLength","ReadLength"])

	sys.stderr.write("Tabulating results for each marker... \n")
	for strMarker in dictHitsForMarker.keys():
		iHits = dictHitsForMarker[strMarker]
		iMarkerNucs = dictMarkerLen[strMarker]*3
		if strCentCheck=="Y":
			strProtFamily = strMarker
			dCount = iHits / float(iMarkerNucs)
			iPossibleHitSpace = float(iMarkerNucs)
		else:
			# Correction factor, since we only require dAlnLength of the reads to align. This results in (1-p)*2 Extra Read len on each side
			dPctAdditionalTargetSeq = ((1.0 - dAlnLength)*2.0)*dReadLength

			#Possible Hit Space = Think of this as the "effective" length of our marker.
			#Any reads starting in the possible hit space should get a valid hit on the marker, anything else will not.
			# Hits/Possible Hit Space ~ Hits/Nucleotide. We can't do just nucleotides because we have special rules for matching markers,
			# depending on their length.

			if (iMarkerNucs > (dReadLength*dAlnLength)):
				iPossibleHitSpace = iMarkerNucs + dPctAdditionalTargetSeq -(dReadLength-1)
				#iPossibleHitSpace = [start-dPctAdditionalTargetSeq,end- (trusted read length -1)]
				#  Any reads dPctAdditionalTargetSeq to the left of the marker will have just enough of the marker to be valid hit.
				#  Same for anything after that, until you get to (end-(readlength-1). Anything after that has some overlap, but not enough to be a valid hit.
			else:
				iPossibleHitSpace = dReadLength-iMarkerNucs -1
				#iPossibleHitSpace = [start-(dReadlength-iMarkerNucs),end- (trusted read length -1)]
				#  Any reads after (start-(dReadlength-iMarkerNucs)) will completely overlap the marker, and be a valid hit.
				#  Same for anything after that, until you get to (start). Anything after that has some overlap, but not enough to be a valid hit.
				#  Try this out as  [iMarkerNucs + 2*(dReadLength-iMarkerNucs) -(dReadLength-1)]
				#
				# Just in case, this one worked well.... iPossibleHitSpace = iMarkerNucs + 2*(dReadLength-iMarkerNucs) -(dReadLength-1)

			dCount = iHits/float(iPossibleHitSpace)

			mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
			strProtFamily = mtchProtStub.group(1)

		if iWGSReads >0:
			dCount =  dCount *  (1e12 / iWGSReads)
		else:
			dCount = 0
			sys.stderr.write("WARNING: 0 Reads found in file:" + strFile )
		tupCount = (strProtFamily,strMarker, dCount,dictHitsForMarker[strMarker],dictMarkerLen[strMarker],dReadLength,iPossibleHitSpace)
		atupMarkerCounts.append(tupCount)

	ProcessHitData(atupMarkerCounts, strMarkerResults=strMarkerResults,strFamFile = strResults)

	return
