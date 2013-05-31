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

c_strUSEARCH	= "usearch"

def MakedbUSEARCH ( strMarkers, strDBName):
	p = subprocess.check_call([c_strUSEARCH, "--makeudb_usearch", strMarkers,
	"--output", strDBName])

	return

def RunUSEARCH ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID):

	strFields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql"

	subprocess.check_call(["time","-o", strMarkers + ".time",
		c_strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strBlastOut,"--userfields", strFields,
		"--threads", str(iThreads)])



def StoreHitCounts(strBlastOut,strValidHits,dictHitsForMarker,dictMarkerLen,dictHitCounts,dID,strCentCheck):
	#strBlastOut - BLAST-formatted output from USEARCH
	#strValidHits - File of BLAST hits that meet ShortBRED's ID and Length criteria. Mainly used for evaluation/debugging.
	#dictMarkerLen - Contains each marker/centroid length
	#dictHitCounts - Contains each family's hit count

	csvwHits = csv.writer( open(strValidHits,'a'), csv.excel_tab )

	#Go through the usearch output, for each prot family, record the number of valid hits
	for aLine in csv.reader( open(strBlastOut), csv.excel_tab ):


		iAln = min(dictMarkerLen[aLine[1]] ,args.iAlnMax)

		#If using centroids (Typically only used for evaluation purposes.)....
		if strCentCheck=="Y":
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
	return


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

	return
