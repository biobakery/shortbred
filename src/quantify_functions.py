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
# (contact Jim Kaminski, jjk451@mail.harvard.edu, Jingjing Tang, jatangne@gmail.com).
#####################################################################################

"""
ToDo:
    * Add a function to count tblastn hits
"""


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
    strVersion = strOutput.decode('utf-8').strip().split(" ")[1].split("_")[0]
    return strVersion
    
def CompareVersions(strVersion1,strVersion2):
    """ This function compares two versions to see which is newer.
    It expects the format 'vXX.XX.XX' The number of X's can vary.  
    It returns 1 if V1 is newer, -1 if V2 is newer, 0 if they are the same."""
    
    aiV1,aiV2 = [list(map(int,x.replace("v","").split("."))) for x in [strVersion1,strVersion2]]
    #print(str(aiV1),str(aiV2))
    return ((aiV1 > aiV2) - (aiV1 < aiV2))

def MakedbUSEARCH ( strMarkers, strDBName,strUSEARCH):
    # This functions calls usearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strUSEARCH, "--makeudb_usearch", strMarkers, "--output", strDBName])
	return

def MakedbRapsearch2 ( strMarkers, strDBName,strPrerapPath):
    # This functions calls prerapsearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strPrerapPath, "-d", strMarkers,"-n", strDBName])
	return

def MakedbBLASTnuc ( strMarkers, strDBName,strGenome,dirTmp):
	print("Making blastdb in " + dirTmp)
	# This functions calls usearch to make a database of the ShortBRED markers.
	p = subprocess.check_call([strMarkers, "-in", strGenome, "-out", strDBName,
		"-dbtype", "nucl", "-logfile", dirTmp + os.sep + "blast_nuc_db.log"])
	return

def MakedbDIAMOND (strMarkers, strDBName, strDIAMOND):
    # This function calls diamond to make a database of the ShortBRED markers.
    p = subprocess.check_call([strDIAMOND, "makedb", "--in", strMarkers, "--db", strDBName])
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


def RunUSEARCH (strWGS,strSearchOut , strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):
    # Calls usearch, strFields specifies the output format from usearch.

	strFields = "query+target+id+alnlen+ql+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql"

	subprocess.check_call([
#		"time","-o", str(dirTmp) + os.sep + os.path.basename(strWGS) + ".time",
		strUSEARCH, "--usearch_local", strWGS, "--db", strDB,
		"--id", str(dID),"--userout", strSearchOut ,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])
  
def RunRAPSEARCH2 (strWGS,strSearchOut , strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strRAPSEARCH2):
    # Calls rapsearch2. It currently does not output the length of the wgs reads, which we use use to 1) filter low reads and .

    with open(strWGS,"r") as streamSeq:
        p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strSearchOut ],stdin=streamSeq,stdout=PIPE)
        p.communicate()
        return
    

    """
    with open(strSearchOut ,'w') as fileBlastOut:
        #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-v",str(iAccepts),"-z", str(iThreads),"-b", "0","-u","1"],stdin=PIPE,stdout=PIPE)
        sys.stderr.write("Running rapsearch2, writing output to " + strSearchOut  + "\n")
        with open(strWGS,"r") as streamSeq:
            p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strSearchOut ],stdin=PIPE,stdout=fileBlastOut)
            for seq in SeqIO.parse(streamSeq, "fasta")	:        
                #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-v",str(iAccepts),"-z", str(iThreads),"-b", "0","-u","1"],stdin=streamSeq,stdout=PIPE)
                #p = Popen([strRAPSEARCH2,"-q","stdin","-d", strDB,"-o",strSearchOut ],stdin=streamSeq,stdout=PIPE)
 
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
	


def RunUSEARCHGenome ( strGenome, strSearchOut , strDB,iThreads,dID, dirTmp, iAccepts, iRejects,strUSEARCH):

	strFields ="query+target+id+alnlen+ql+mism+opens+qlo+qhi+tlo+thi+evalue+bits+ql"

	subprocess.check_call([
#		"time","-o", str(dirTmp) + os.sep + os.path.basename(strWGS) + ".time",
		strUSEARCH, "--usearch_local", strGenome, "--db", strDB,
		"--id", str(dID),"--userout", strSearchOut ,"--userfields", strFields,"--maxaccepts",str(iAccepts),
		"--maxrejects",str(iRejects),"--threads", str(iThreads)])

def RunTBLASTN ( strTBLASTN, strDB,strGenome, strSearchOut , iThreads):



	astrBlastParams = ["-outfmt", "6 sseqid qseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen",
	"-matrix", "BLOSUM62", "-ungapped",
		"-comp_based_stats","F","-window_size","0",
		"-xdrop_ungap","1","-evalue","1e-3",
		"-max_target_seqs", "1000000",
		"-num_threads",str(iThreads)]


	subprocess.check_call(
		[strTBLASTN, "-db", strDB,"-query", strGenome,"-out",strSearchOut ] +  astrBlastParams)

def RunDIAMONDx (strDIAMOND, strDB, strWGS, strDiamondOut, iThreads):
    astrDIAMONDParams = ["--outfmt", "6","qseqid","sseqid","pident","length","qlen","mismatch","gaps","qstart","qend","sstart","send",
                         "evalue","bitscore","qlen", "--matrix", "BLOSUM62",
                       "--comp-based-stats","0","--window","0","--ungapped-score","20",
                       "--evalue","1e-8","--rank-ratio","0","--shapes","0",
                       "--max-target-seqs", "0", 
                       "--threads",str(iThreads)]
    # "--shape-mask","11111" this param is removed temporarily, otherwise, some of self-hits would disappear
    
    cmdRun = [
#        "time", "-o", dirTime + os.sep +"goisearch.time",
        strDIAMOND,"blastx", "--query", strWGS, "--db", strDB,
        "--out", strDiamondOut] + astrDIAMONDParams
            
    #sys.stderr.write(" ".join(cmdRun))
    subprocess.check_call(cmdRun)    
    return

def RunDIAMONDp (strDIAMOND, strDB, strGenome, strDiamondOut, iThreads):
    astrDIAMONDParams = ["--outfmt", "6","qseqid","sseqid","pident","length","mismatch","gaps","qstart","qend","sstart","send",
                         "evalue","bitscore","qlen", "--matrix", "BLOSUM62",
                       "--comp-based-stats","0","--window","0","--ungapped-score","20",
                       "--evalue","1e-8","--rank-ratio","0","--shapes","0",
                       "--max-target-seqs", "0",
                       "--threads",str(iThreads)]
    # "--shape-mask","11111" this param is removed temporarily, otherwise, some of self-hits would disappear
    
    cmdRun = [
#        "time", "-o", dirTime + os.sep +"goisearch.time",
        strDIAMOND,"blastp", "--query", strWGS, "--db", strDB,
        "--out", strDiamondOut] + astrDIAMONDParams
            
    #sys.stderr.write(" ".join(cmdRun))
    subprocess.check_call(cmdRun)    
    return

def Median(adValues):
	adValues.sort()
	iLen = len(adValues)
	if iLen % 2==0:
		i = int(iLen/2)
		dMedian = float(adValues[i] + adValues[i-1])/2.0
	else:
		dMedian = adValues[int(math.floor(iLen/2))]
	return dMedian

def StoreHitCounts(strSearchOut ,strValidHits,dictHitsForMarker,dictMarkerLen,dictCountsForMarker, dID,strCentCheck,dAlnLength,iMinReadBP,iAvgMarkerAA,strSearchMethod, strShortBREDMode="wgs",iAlnCentroids=30,version_control = 1):
# Reads in the USEARCH output (strSearchOut ), marks which hits are valid (id>=dID &
# len >= min(95% of read,dictMarkerLen[Marker]) and adds to count in dictHitsForMarker[strMarker].
# Valid hits are also copied to the file in strValidHits. strCentCheck is used to flag centroids,
# and handle their counting
    
    

    if strSearchMethod == "usearch":
        len_control = CompareVersions(version_control,c_vstrUsearchForAAValue)
    else:
        len_control = version_control
    
    # If the version in use is newer than 6.0.307,     
    with open(strValidHits, 'a') as csvfileHits:
        csvwHits = csv.writer( csvfileHits, csv.excel_tab )
        sys.stderr.write("Processing %s results... \n" % strSearchMethod.upper())
        #Go through the usearch output, for each prot family, record the number of valid
        with open(strSearchOut , 'r') as csvfileSearch:
            for aLine in csv.reader( csvfileSearch, delimiter='\t' ):
                strMarker 	= aLine[1]
                dHitID		= aLine[2]      # percentage of alignment
                iAlnLen     = int(aLine[3])   # alignment length
                RLen        = int(aLine[4])          # read length
                if RLen == 0:
                    continue
#                iReadLenAA = int(aLine[8]) - int(aLine[7]) + 1
                
                # USearch versions past 6.0.307 report query length in AA space, 6.0.307 and prior versions report in nuc space.
                # This is not an issue in AA to AA comparisons, so seraching "annnotated_genomes" is fine. 
#                if strShortBREDMode!="wgs" or len_control != 1:
#                    iReadLenAA = iReadLenAA 
#                else:
#                    iReadLenAA = int(round(iReadLenAA/3))                    
                
                # A valid match must be as long as 95% of the read or the full marker. Should be changed'
                # (Note that this in AA's.)
                #If using centroids (Typically only used for evaluation purposes.)....
                if strCentCheck=="Y":
                    strProtFamily = strMarker
                    
                    if ( (int(iAlnLen)>= iAlnCentroids) and ( float(dHitID)/100) >= dID):
                        dictHitsForMarker[strProtFamily] = dictHitsForMarker.setdefault(strProtFamily,0) + 1
                        csvwHits.writerow( aLine )
                    
                
                #If using ShortBRED Markers (and not centroids)...
                else:
                    iAlnMin = min(dictMarkerLen[strMarker],math.floor(RLen*dAlnLength/3)) # marker length
                    #Get the Family Name
                    mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
                    strProtFamily = mtchProtStub.group(1)
                    """
                    if(float(dHitID)>.9):
                        sys.stderr.write( " ".join([str(iAlnLen),str(iAlnMin),str(dHitID),str(dID),str(iReadLenAA),str(iMinMarkerAA)] ) +"\n")
                    """
                    #If hit satisfies criteria, add it to counts, write out data to Hits file
                    if (int(iAlnLen)>= iAlnMin and (RLen>=iMinReadBP) and (float(dHitID)/100) >= dID):
                        # criteria (iReadLenAA >= iMinMarkerAA) is changed to RLen > args.iMinReadBP (both of them are read lengths)
                        # iReadLenAA is the length of read alignment, is it right to compare it with 
                        # iMinmarkerAA = =int(math.floor(args.iMinReadBP/3)) where iMinReadBP is the lower bound for read lengths 
                        # that shortbred will process
                        
                        #Add 1 to count of hits for that marker, and family
                        iMarkerNucs = dictMarkerLen[strMarker]*3
                        dPctAdditionalTargetSeq = ((1.0 - dAlnLength)*2.0)*RLen
                        
                        if (iMarkerNucs > (RLen*dAlnLength)):
                            iPossibleHitSpace = iMarkerNucs + dPctAdditionalTargetSeq -RLen+1
                        
                        else:
                            iPossibleHitSpace = RLen-iMarkerNucs -1

                        dictHitsForMarker[aLine[1]] = dictHitsForMarker.setdefault(aLine[1],0) + 1
                        dictCountsForMarker[aLine[1]] += float(1/iPossibleHitSpace)
                        csvwHits.writerow( aLine )
    return



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

def CalculateCounts(strResults,strMarkerResults, dictHitsForMarker, dictCountsForMarker,dictMarkerLen,iWGSReads,strCentCheck,
dAlnLength,strFile):
    #strResults - Name of text file with final ShortBRED Counts
    #strSearchOut - BLAST-formatted output from USEARCH/RAPSEARCH2/DIAMOND
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
            dCount = iHits/(float(iPossibleHitSpace)/1000)
        else:
            dCount = dictCountsForMarker[strMarker]
            
        mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',strMarker)
        strProtFamily = mtchProtStub.group(1)
        
        if iWGSReads >0:
            dCount =  dCount /  (iWGSReads / 1e9 )
        
        else:
            dCount = 0
            sys.stderr.write(str(strFile))            
            sys.stderr.write("WARNING: 0 Reads found in file:" + strFile )
        
        tupCount = (strProtFamily,strMarker, dCount,dictHitsForMarker[strMarker],dictMarkerLen[strMarker])
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


