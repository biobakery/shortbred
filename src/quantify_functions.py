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
from subprocess import Popen, PIPE,STDOUT
import csv
import re
import sys
import math
import os
import io

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

#c_iAlnCentroids = 30
# This used to be a constant, but we have changed it to allow for different
# values when working with centroids.

c_vstrUsearchForAAValue = "v6.0.307"
# Explained in CheckUSEARCH function.

def CheckUSEARCH(strUSEARCH):
    """ This function returns the version of usearch. 
    In searches with a nucleotide query and amino acid db, versions v6.1.544
    and beyond report query length in nucleotides, and versions v6.0.307 and
    older report query length in amino acids. 
    
    We pass the version to StoreHitCounts, so that the correct amino acid
    length is used. """
    strOutput = subprocess.check_output([strUSEARCH, "--version"])
    strVersion = strOutput.strip().split(" ")[1].split("_")[0]
    return strVersion
    
def CompareVersions(strVersion1,strVersion2):
    """ This function compares two versions to see which is newer.
    It expects the format 'vXX.XX.XX' The number of X's can vary.  
    It returns 1 if V1 is newer, -1 if V2 is newer, 0 if they are the same."""
    
    aiV1,aiV2 = [map(int,x.replace("v","").split(".")) for x in [strVersion1,strVersion2]]
    #print(str(aiV1),str(aiV2))
    return cmp(aiV1,aiV2)

def MakedbUSEARCH ( strMarkers, strDBName,strUSEARCH):
    # This functions calls usearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strUSEARCH, "--makeudb_usearch", strMarkers,"--output", strDBName])
	return

def MakedbRapsearch2 ( strMarkers, strDBName,strPrerapPath):
    # This functions calls prerapsearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strPrerapPath, "-d", strMarkers,"-n", strDBName])
	return

def MakedbBLASTnuc ( strMakeBlastDB, strDBName,strGenome,dirTmp):
	print "Making blastdb in " + dirTmp
	# This functions calls usearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strMakeBlastDB, "-in", strGenome, "-out", strDBName,
		"-dbtype", "nucl", "-logfile", dirTmp + os.sep + "blast_nuc_db.log"])
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

def MakeDictFamilyCounts (strMarkers,strFamilyOut):
	# Counts up the number of markers each protein family has,
	# saves it to dictFamMarkerCounts.

	dictFamMarkerCounts = {}
	sys.stderr.write("Calculating markers per family... \n")
	for seq in SeqIO.parse(strMarkers, "fasta"):
		mtchFam = re.search(r'^(.*)_[TJQ]M.*',seq.id)
		if(mtchFam):
			strFam = str(mtchFam.group(1)).strip()
			if strFam in dictFamMarkerCounts:
				dictFamMarkerCounts[strFam] = dictFamMarkerCounts[strFam]+1
			else:
				dictFamMarkerCounts[strFam] = 1
	return dictFamMarkerCounts

def CalcFinalCount (dictORFMatches,dictFamMarkerCounts,bUnannotated,dPctORFScoreThresh,dPctMarkerThresh):
    # Takes two dictionaries, each have protein families has the keys.
	# One has the number of markers hitting the ORF, the other has all possible markers.

    # Make two paramaters her avaialble as arguments:
    # 1) dPctORFScoreThresh, which controls how much a share of an ORF score a fam must receive
    # 2) dPctMarkerThresh, which controls the threshhold for even being counted as a hit
	

	aaCounts = []
	aaFinalCounts = []


	for strFam in dictORFMatches:
		dScore = dictORFMatches[strFam] / float(dictFamMarkerCounts[strFam])
		#print strFam,dScore,dPctMarkerThresh
		if (dScore < dPctMarkerThresh):
			dScore = 0
			#print strFam,dScore,dPctMarkerThresh
  
		aFamScore = [strFam,dScore]
		aaCounts.append(aFamScore)

	dSum = sum(list(zip(*aaCounts))[1])

	# For annotated genomes, we throw out scores that comprise less than (dThresh)
	# percent of the total for the ORF
	if (bUnannotated==False):
		aaAboveThresh = []
		for aFamScore in aaCounts:
			if dSum >0:
				if (aFamScore[1] / dSum)>=dPctORFScoreThresh:
					aaAboveThresh.append(aFamScore)
					#print aaAboveThresh, dSum

		for aFamScore in aaAboveThresh:
			aNewScore = [aFamScore[0],aFamScore[1] * (aFamScore[1]/dSum) ]
			aaFinalCounts.append(aNewScore)

	# For unannotated genomes, multiple families can hit to a contig, and we
	# want to count all of them. We do not apply the threshold.
	else:
		aaFinalCounts = aaCounts
	return aaFinalCounts

	"""
	Example:

	Before
    	FamA	0.85
		FamB	0.32

	After
		FamA	0.62
		FamB	0.09

	"""

def NormalizeGenomeCounts (strValidHits,dictFamCounts,bUnannotated,dPctORFScoreThresh=.1,dPctMarkerThresh=.1):
	dictFinalCounts = {}
	for strFam in dictFamCounts.keys():
		dictFinalCounts[strFam] = 0
	dictORFMatches = {}

	# Make dictionart where
	# strORF ~ set (Marker1, Marker2, ...)
	with open(strValidHits, 'r') as csvfileHits:
		for aLine in csv.reader( csvfileHits, delimiter='\t' ):

			strORF = aLine[0]
			strMarker = aLine[1]


			if strORF in dictORFMatches:
				dictORFMatches[strORF] = dictORFMatches[strORF] + [strMarker]

			else:
				dictORFMatches[strORF] = [strMarker]


	for strORF in sorted(dictORFMatches.keys()):
		if(bUnannotated==False):
			astrMatches = set(dictORFMatches[strORF])
		else:
			astrMatches = dictORFMatches[strORF]

		dictFamMatches = {}


		# get count of matches to each family
		for strFam in astrMatches:
			mtchFam = re.search(r'^(.*)_[TJQ]M.*',strFam)
			if(mtchFam):
				strFam = str(mtchFam.group(1)).strip()

				if strFam in dictFamMatches:
					dictFamMatches[strFam] = dictFamMatches[strFam]+1
				else:
					dictFamMatches[strFam] = 1

		# Normalize Counts
		aaCount = CalcFinalCount (dictFamMatches,dictFamCounts,bUnannotated,dPctORFScoreThresh,dPctMarkerThresh)

		for aFamScore in aaCount:
			dictFinalCounts[aFamScore[0]] = dictFinalCounts[aFamScore[0]] + aFamScore[1]

	return dictFinalCounts


def RunUSEARCH ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):
    # Calls usearch, strFields specifies the output format from usearch.

	strFields = "query+target+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts"

	subprocess.check_call([
#		"time","-o", str(dirTmp) + os.sep + os.path.basename(strMarkers) + ".time",
		strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strBlastOut,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])
  
def RunRAPSEARCH2 ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strRAPSEARCH2):
    # Calls rapsearch2. It currently does not output the length of the wgs reads, which we use use to 1) filter low reads and .

    with open(strWGS,"r") as streamSeq:
        p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strBlastOut],stdin=streamSeq,stdout=PIPE)
        p.communicate()
    	return
    

    """
    with open(strBlastOut,'w') as fileBlastOut:
        #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-v",str(iAccepts),"-z", str(iThreads),"-b", "0","-u","1"],stdin=PIPE,stdout=PIPE)
        sys.stderr.write("Running rapsearch2, writing output to " + strBlastOut + "\n")
        with open(strWGS,"r") as streamSeq:
            p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strBlastOut],stdin=PIPE,stdout=fileBlastOut)
            for seq in SeqIO.parse(streamSeq, "fasta")	:        
                #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-v",str(iAccepts),"-z", str(iThreads),"-b", "0","-u","1"],stdin=streamSeq,stdout=PIPE)
                #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strBlastOut],stdin=streamSeq,stdout=PIPE)
 
                 print seq.format("fasta")
 
                 p.communicate(seq.format("fasta"))
                 #seq = next(streamSeq)
                 #print(seq.format("fasta"))
                 #strRapOut = p.communicate(input=seq.format("fasta"))
                 #p.stdin.write(seq.format("fasta"))
                 #p.stdin.close()
                 #strRapOut = p.stdout.read()
                 #print(strRapOut)
    """
	


def RunUSEARCHGenome ( strMarkers, strWGS,strBlastOut, strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):

	strFields = "target+query+id+alnlen+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql+tl+qs+ts"

	subprocess.check_call([
#		"time","-o", str(dirTmp) + os.sep + os.path.basename(strMarkers) + ".time",
		strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strBlastOut,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])

def RunTBLASTN ( strTBLASTN, strDB,strMarkers, strBlastOut, iThreads):



	astrBlastParams = ["-outfmt", "6 sseqid qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
	"-matrix", "PAM30", "-ungapped",
		"-comp_based_stats","F","-window_size","0",
		"-xdrop_ungap","1","-evalue","1e-3",
		"-max_target_seqs", "1000000",
		"-num_threads",str(iThreads)]


	subprocess.check_call(
		[strTBLASTN, "-db", strDB,"-query", strMarkers,"-out",strBlastOut] +  astrBlastParams)

def Median(adValues):
	adValues.sort()
	iLen = len(adValues)
	if iLen % 2==0:
		i = int(iLen/2)
		dMedian = float(adValues[i] + adValues[i-1])/2.0
	else:
		dMedian = adValues[int(math.floor(iLen/2))]
	return dMedian

def StoreHitCountsRapsearch2(strBlastOut,strValidHits,dictHitsForMarker,dictMarkerLen,dictHitCounts,dID,strCentCheck,dAlnLength,iMinReadAA,iAvgReadAA,iAlnCentroids=30,strUSearchOut=False):
# Reads in the RAPSEARCH2 output (strBlastOut), marks which hits are valid (id>=dID &
# len >= min(95% of read,dictMarkerLen[Marker]) and adds to count in dictHitsForMarker[strMarker].
# Valid hits are also copied to the file in strValidHits. strCentCheck is used to flag centroids,
# and handle their counting

# Some small changes are made for Rapsearch2 output in this function.
# *add ".m8", to the filename 
# *skip the first 5 lines
#* rearrange columns.

    strBlastOut = strBlastOut + ".m8"
    iSkip = 5

    with open(strValidHits, 'a') as csvfileHits:
        csvwHits = csv.writer( csvfileHits, csv.excel_tab )
        sys.stderr.write("Processing RAPSEARCH2 results... \n")
        #Go through the Rapsearch2 output, for each prot family, record the number of valid hits

        with open(strBlastOut, 'r') as csvfileBlast:
            csvReader = csv.reader( csvfileBlast, delimiter='\t' )
            for i in range(iSkip):
    			next(csvReader)
            for aLine in csvReader:
                strMarker       = aLine[1]
                dHitID          = aLine[2]
                iAlnLen     = int(aLine[3])
                if (strUSearchOut):
                    iReadLenAA  = int(aLine[12])
                else:
                    iReadLenAA = abs(int(aLine[7]) - int(aLine[6]))/3.0
    
                # A valid match must be as long as 95% of the read or the full marker.
                # (Note that this in AA's.)
    
                #If using centroids (Typically only used for evaluation purposes.)....
                if strCentCheck=="Y":
                    strProtFamily = strMarker
                    
                    if ( (int(iAlnLen)>= iAlnCentroids) and ( float(dHitID)/100) >= dID):
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
                    # NOTE: Rapsearch2 does not currently output the original query read length, so we do *NOT* filter based on read length.
                    if (int(iAlnLen)>= iAlnMin and (float(dHitID)/100) >= dID):
    
                        #Add 1 to count of hits for that marker, and family
                        dictHitsForMarker[aLine[1]] = dictHitsForMarker.setdefault(aLine[1],0) + 1
                        dictHitCounts[strProtFamily] = dictHitCounts.setdefault(strProtFamily,0) +1
    
                        csvwHits.writerow( aLine )
        return


def StoreHitCounts(strBlastOut,strValidHits,dictHitsForMarker,dictMarkerLen,dictHitCounts,dID,strCentCheck,dAlnLength,iMinReadAA,iAvgReadAA,strVersionUSEARCH,strShortBREDMode="wgs",iAlnCentroids=30,strUSearchOut=True):
# Reads in the USEARCH output (strBlastOut), marks which hits are valid (id>=dID &
# len >= min(95% of read,dictMarkerLen[Marker]) and adds to count in dictHitsForMarker[strMarker].
# Valid hits are also copied to the file in strValidHits. strCentCheck is used to flag centroids,
# and handle their counting
    
	# If the version in use is newer than 6.0.307,     
	CompareVersions(strVersionUSEARCH,c_vstrUsearchForAAValue)


	with open(strValidHits, 'a') as csvfileHits:
		csvwHits = csv.writer( csvfileHits, csv.excel_tab )
		sys.stderr.write("Processing USEARCH results... \n")
		#Go through the usearch output, for each prot family, record the number of valid

		with open(strBlastOut, 'r') as csvfileBlast:
			for aLine in csv.reader( csvfileBlast, delimiter='\t' ):
				strMarker 	= aLine[1]
				dHitID		= aLine[2]
				iAlnLen     = int(aLine[3])
				if (strUSearchOut):
					iReadLenAA  = int(round(int(aLine[12])))
				else:
					iReadLenAA = int(aLine[7]) - int(aLine[6])
				
				# USearch versions past 6.0.307 report query length in AA space, 6.0.307 and prior versions report in nuc space.
				# This is not an issue in AA to AA comparisons, so seraching "annnotated_genomes" is fine. 
 				if (strShortBREDMode!="wgs" or CompareVersions(strVersionUSEARCH,c_vstrUsearchForAAValue) != 1):    
					iReadLenAA = iReadLenAA 
 				else:
					iReadLenAA = int(round(iReadLenAA/3))                    
                    

				# A valid match must be as long as 95% of the read or the full marker.
		        # (Note that this in AA's.)
				#If using centroids (Typically only used for evaluation purposes.)....
				if strCentCheck=="Y":
					strProtFamily = strMarker

					if ( (int(iAlnLen)>= iAlnCentroids) and ( float(dHitID)/100) >= dID):
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
	atupZipped = list(zip(*atupCurFamData))

	# Family Stats
	try:
		dMedian = Median(list(atupZipped[2]))
		iHits = sum(atupZipped[3])
		iMarkerLength = sum(atupZipped[4])
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
		iHits = dictHitsForMarker.get(strMarker,0)
		iMarkerNucs = dictMarkerLen[strMarker]*3
		if strCentCheck=="Y":
			strProtFamily = strMarker
			dCount = iHits / (float(iMarkerNucs)/1000)
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

			dCount = iHits/(float(iPossibleHitSpace)/1000)

			mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
			strProtFamily = mtchProtStub.group(1)

		if iWGSReads >0:
			dCount =  dCount /  (iWGSReads / 1e6 )
		else:
			dCount = 0
			sys.stderr.write("WARNING: 0 Reads found in file:" + strFile )
		tupCount = (strProtFamily,strMarker, dCount,dictHitsForMarker[strMarker],dictMarkerLen[strMarker],dReadLength,iPossibleHitSpace)
		atupMarkerCounts.append(tupCount)

	ProcessHitData(atupMarkerCounts, strMarkerResults=strMarkerResults,strFamFile = strResults)

	return atupMarkerCounts

####################################################################################################
def BayesUpdate(atupCounts,strBayesResults,strBayesLog,astrQMs,dictQMPossibleOverlap,dictType):
	# This function uses Bayesian updating to adjust the final ShortBRED counts
	# for quasi markers. We consider true markers and junction markers to provide
	# good estimates of the abundance of proteins of interest. By subtracting these
	# values from the estimated QM counts, we obtain better estimates of their
	# true values.


	strCurFam = ""
	atupCurFamData = []
	dictFamCounts = {}
	dictMarkerCounts = {}


	for tupRow in atupCounts:

		strFam = tupRow[0]
		dictMarkerCounts[tupRow[1]] = tupRow[2]

		# For each family, get Median ShortBRED count.
		if strFam != strCurFam and strCurFam!="":
			# If start a new family, print out the median for the last one.
			atupZip = zip(*atupCurFamData)
			dMed = Median(list(atupZip[2]))
			dictFamCounts[strCurFam] = dMed
			atupCurFamData = [tupRow]
			strCurFam = strFam
		elif strCurFam=="":
			# For the first fam, initialize the family and current data.
			strCurFam=strFam
			atupCurFamData=[tupRow]
		else:
			# If still on same family, add row to current data.
			atupCurFamData+=[tupRow]

		# Add count for the last family.
		atupZip = zip(*atupCurFamData)
		dMed = Median(list(atupZip[2]))
		dictFamCounts[strCurFam] = dMed

	# The BayesLog was used for debugging purposes, this code also updates
	# the value for the quasi-markers.

	with open(strBayesLog, "a") as BayesLog:
		BayesLog.write("Total Families:" +str(len(dictFamCounts)) + "\n")

		for strMarker in astrQMs:
			iQM = 0
			iJM = 0
			iTM = 0
			dSubtract =0
			BayesLog.write(strMarker+": ")
			for strFam in dictQMPossibleOverlap[strMarker]:
				BayesLog.write(strFam + " " + str(dictFamCounts[strFam]) + "\n")
				if dictType[strFam]=="TM":
					iTM+=1
					dSubtract+=dictFamCounts[strFam]
				elif dictType[strFam]=="JM":
					iJM+=1
					dSubtract+=dictFamCounts[strFam]
				elif dictType[strFam]=="QM":
					strQMFam = strFam
					iQM+=1

			BayesLog.write("\n")
			astrValues=[str(iTM),str(iJM),str(iQM)]
			dTotal = dictMarkerCounts[strMarker] - dSubtract
			if dTotal<0:
				dTotal = 0
			BayesLog.write("\t".join(astrValues) + "\n")
			BayesLog.write("Final Value: " +str(dTotal) + "\n")
			dictFamCounts[strQMFam] = dTotal

	with open(strBayesResults, "w") as fOut:
		fOut.write("Family" + "\t" + "Count" + "\n")
		for strFam in dictFamCounts.keys():
			fOut.write(strFam + "\t" + str(dictFamCounts[strFam]) + "\n")


	return


