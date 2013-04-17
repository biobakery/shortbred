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

import re
import sys
import csv
import argparse
import os
import subprocess
import glob
import math
import shutil

import Bio
from Bio import AlignIO
from Bio import Align
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

###############################################################################
def GetCDHitMap ( fileCluster, strTxtMap):
    setGenes = set()
    strFam = ""
    iLine = 0

    dictFams = {}

    for astrLine in csv.reader(open(fileCluster),delimiter='\t'):
		iLine+=1
		mtch = re.search(r'^>',astrLine[0])
		if (mtch):
		    if (len(setGenes)>0 and strFam!=""):
				for strGene in setGenes:
				    dictFams[strGene] = strFam
				setGenes=set()
				strFam = ""
		else:
		    mtchGene = re.search(r'>(.*)\.\.\.',str(astrLine[1]))
		    strGene = mtchGene.groups(1)[0]
		    setGenes.add(strGene.strip())
		    mtchFam = re.search(r'\.\.\.\s\*',str(astrLine[1]))
		    if mtchFam:
				strFam = strGene.strip()

    for strGene in setGenes:
				    dictFams[strGene.strip()] = strFam.strip()


    f = open(strTxtMap, 'w')
    for prot, fam in sorted(dictFams.items(), key = lambda(prot, fam): (fam,prot)):
		f.write(fam + "\t" + prot + "\n")
    f.close()

###############################################################################
def MakeFamilyFastaFiles ( strMapFile, fileFasta, dirOut):
#Makes a fasta file containing genes for each family in dictFams.
    dictFams = {}
    for strLine in csv.reader(open(strMapFile),delimiter='\t'):
		dictFams[strLine[1]]=strLine[0]


    dictgeneFamilies = {}

    for gene in SeqIO.parse(fileFasta, "fasta"):
		strFamily = dictFams.get(gene.id,"unassigned")
		aGene = dictgeneFamilies.get(strFamily,list())
		aGene.append(gene)
		dictgeneFamilies[strFamily] = aGene

    for key in dictgeneFamilies.keys():
		strFamFile = dirOut + os.sep + key +".faa"
		f = open(strFamFile, 'w')
		SeqIO.write(dictgeneFamilies[key], strFamFile, "fasta")
		f.close()

###############################################################################
def ClusterFams(dirClust, dCLustID, strOutputFile, dThresh):
#Clusters all of the family files made by MakeFamilyFastaFiles.

    dirFams = dirClust + os.sep + "fams"
    dirCentroids = dirFams+os.sep+"centroids"
    dirUC = dirFams+ os.sep + "uc"

    if not os.path.exists(dirClust):
		os.makedirs(dirClust)
    if not os.path.exists(dirFams):
		os.makedirs(dirFams)
    if not os.path.exists(dirCentroids):
		os.makedirs(dirCentroids)

    if not os.path.exists(dirUC):
		os.makedirs(dirUC)


    print dirCentroids
    print glob.glob(dirFams+os.sep+'*.faa')
    for fileFasta in glob.glob(dirFams+os.sep+'*.faa'):
		print "The file is ", fileFasta
		fileClust = dirCentroids + os.sep + os.path.basename(fileFasta)
		fileAlign = dirFams + os.sep + os.path.basename(fileFasta)+".aln"
		strSeqID = os.path.basename(fileFasta)
		strSeqID = strSeqID.replace(".faa", "")

		iSeqCount = 0
		#Count seqs, if more than one, then align them
		for seq in SeqIO.parse(fileFasta, "fasta"):
		    iSeqCount+=1


		if iSeqCount>1:
		    #Call muscle to produce an alignment
			subprocess.check_call(["muscle", "-in", str(fileFasta), "-out", str(fileAlign)])

            # Use BioPython's "dumb consensus" feature to get consensus sequence
			algnFasta = AlignIO.read(str(fileAlign), "fasta")

			seqConsensus =str(AlignInfo.SummaryInfo(algnFasta).dumb_consensus(threshold=dThresh, ambiguous='X'))
			seqConsensus = SeqRecord(Seq(seqConsensus),id=strSeqID)

			SeqIO.write(seqConsensus, str(fileClust), "fasta")


			"""
			# We previously used EMBOSS-CONS to produce consensus sequences
			# Call cons or em_cons from the EMBOSS package to produce a consensus sequence
			subprocess.check_call(["cons", "-seq", str(fileAlign), "-outseq", str(fileClust)])
			"""
		else:
		    shutil.copyfile(fileFasta,fileClust)



    ageneAllGenes = []

    for fileFasta in glob.glob(dirCentroids+os.sep+'*.faa'):
		for gene in SeqIO.parse(fileFasta, "fasta"):
		    gene.id = os.path.basename(os.path.splitext(fileFasta)[0])
		    ageneAllGenes.append(gene)

    """
    for gene in ageneAllGenes:
		mtch = re.search(r'centroid=(.*)',gene.id)
		if mtch:
		    gene.id = mtch.group(1)
		else:
		    gene.id = os.path.splitext()
    """

    SeqIO.write(ageneAllGenes, strOutputFile, "fasta")

###############################################################################
def printMap(strUCMap,strTxtMap):
#Converts usearch .uc file into two column (FamId, GeneId) file.

    dictGeneMap={}
    for strLine in csv.reader(open(strUCMap),delimiter='\t'):
		if (strLine[0] == "H"):
		    dictGeneMap[strLine[-2]] = strLine[-1]
		elif (strLine[0] == "C"):
		    dictGeneMap[strLine[-2]] = strLine[-2]

    f = open(strTxtMap, 'w')
    for prot, fam in sorted(dictGeneMap.items(), key = lambda(prot, fam): (fam,prot)):
		f.write(fam + "\t" + prot + "\n")
    f.close()

###############################################################################
def getGeneData ( fileFasta):
#Returns dict of form (GeneID, seq)
    dictGeneData = {}

    for gene in SeqIO.parse(fileFasta, "fasta"):
		    dictGeneData.setdefault(gene.id.strip(), str(gene.seq))


    if dictGeneData.has_key(''):
		del dictGeneData['']
    return dictGeneData

##############################################################################
def getOverlapCounts (fileBlast, dIDcutoff, iLengthMin, dLengthMax, iOffset, bSaveHitInfo):

# Makes a dictionary of form (GeneID, [0,0,0,0,1,1...]), where the number indicates
# the number of times an amino acid overlaps with a region in the blast output.

# PARAMETERS
# fileBlast - BLAST-formatted (fmt6) file of hits
# dIDcutoff - Minimum ID for using a hit
# iLengthMin - Lower bound for using a hit in AA's
# dLengthMax - Upper bound for using a hit (AlignmentLength/QueryLength)
# iOffset   - Number of AminoAcids at each end of valid hit which do NOT get +1 to their overlap counts
# bSaveHitInfo - When true, this records the name of the *Target* hit, and the start and end of the hit on the *Query*




#Read in the blast output line by line
#When the program finds a new QueryGene:
    #Add the last gene's aiCounts to dictAAOverlapCounts -- (gene, abWindow)
	#Add the last gene's atupHitInfo to dictOverlapInfo
    #Make a new aiCounts
#When the program finds the same QueryGene:
    #If the region in the blast hit satisfies our length and ID parameters
		#Set the corresponding location in abWindow equal to (count+1)

#Add the last gene when you finish the file

	strCurQuery = ""
	dictAAOverlapCounts = {}
	dictOverlapInfo = {}
	atupHitsInfo = []
	aiCounts =[]
	iGeneCount = 0
	iLine =0


    #sys.stderr.write("The file is " + fileBlast + "\n")
	for aLine in csv.reader( open(fileBlast, 'rb'), delimiter='\t' ):

		iLine+=1

		strQueryID = aLine[0]
		strSubId = aLine[1]
		dIdentity =float(aLine[2])/100
		iAln= int(aLine[3] )
		iMismatch = int(aLine[4])
		iGap = int(aLine[5] )
		iQStart = int(aLine[6] )
		iQEnd = int(aLine[7])
		iSubStart = int(aLine[8] )
		iSubEnd= int(aLine[9] )
		deVal = float(aLine[10] )
		dBit = float(aLine[11])
		iQLength = int(aLine[12])

		#When the current queryID changes, add the last one, and switch to the next gene.
		if strQueryID != strCurQuery:
			if iLine>1:
				dictAAOverlapCounts.setdefault(strCurQuery, aiCounts)
				dictOverlapInfo.setdefault(strCurQuery.strip(), atupHitsInfo)
			iGeneCount = iGeneCount+1
			strCurQuery = strQueryID

			atupHitsInfo = []
			aiCounts = []
			for i in range(iQLength):
				aiCounts.append(0)

		dMatchLength = (iAln) / float(iQLength)


		#If it's valid overlap:
		#   and 1 to each appopriate spot in the array
		#   save the name of the overlapping gene, and its start and end points
		#       Note: The latter region will include AA's outside what is marked as overlap, we still want that info.

		if (dIdentity >= dIDcutoff) and (iAln >=iLengthMin) and (dMatchLength <= dLengthMax) and (strQueryID!=strSubId):
			iOLCount = 0
			for i in range(iQStart-1+iOffset, iQEnd-iOffset):
				aiCounts[i]=aiCounts[i]+1
				iOLCount += 1
			if(bSaveHitInfo == True):
				tupInfo = strSubId, iQStart, iQEnd, iSubStart,iSubEnd, iOLCount
				atupHitsInfo.append(tupInfo)



	dictAAOverlapCounts.setdefault(strCurQuery.strip(), aiCounts)
 	dictOverlapInfo.setdefault(strCurQuery.strip(), atupHitsInfo)
  	iGeneCount = iGeneCount+1

	return dictAAOverlapCounts, dictOverlapInfo
	#tupAAandHits = (dictAAOverlapCounts, dictOverlapInfo)



###########################################################################
def MarkX(dictGenes, dictOverlap):

    for strName in dictGenes:
		for i in range(len(dictGenes[strName])):
		    if dictGenes[strName][i]=="X" or dictGenes[strName][i]=="x":
				dictOverlap[strName][i] = dictOverlap[strName][i] + 999999

    return dictOverlap


###########################################################################
def CheckForMarkers(setGenes, dictKnockOut, iN):
#Take the genes in setNames, look in dictKnockOut (GeneID, [0,0,..])
#to see if they have a region of length N without any overlap.
# Returns set of genes (ID's) with markers.


    fFoundRegion = False
    fHitEnd = False
    iSumCounts = 0
    iMarkers = 0

    setGenesWithMarkers = set()

    for key in setGenes:
		aiWindow = dictKnockOut[key]
		iStart = 0
		while (fFoundRegion == False and fHitEnd == False):
		    if ((iStart+iN) > len(aiWindow)):
				fHitEnd = True
		    iSumCounts = sum(aiWindow[iStart:(iStart+iN)])
		    if (iSumCounts == 0 and fHitEnd==False):
				iMarkers +=1
				fFoundRegion = True
				setGenesWithMarkers.add(key)
		    iStart += 1
		fFoundRegion = False
		fHitEnd = False

    return setGenesWithMarkers
###############################################################################
def CheckForQuasiMarkers(setGenes, dictKnockOut, dictGenes, iN, iThresh, iTotLength):

    #Only run this on the leftover genes
    # "iN" = minimum window length
    #For each one, sum up the values from [0:n], then [1:n+1]...
    #Store these in an array of length (len(gene)-n)
    #Find the minimum value in this array
    #Take its index
    #Your window is [index:index+n]

    #Get the appropriate string from dictGOIgenes
    #add it to dictQM


    #Return dictQM with these windows

    #A QM_tuple has 4 values (name, markers, overlapvalue,iOriginalLength)


    atupQM = []



    for key in setGenes:
		aiWindow = dictKnockOut[key]

		#Take some function of each overlap count, reduce influence of outliers

		import math

		adAdjWindow = [math.pow(x,(1/4.0)) for x in aiWindow]

		iStart = 0
		dMin = 0
		dSumCounts = 0
		fHitEnd = False
		adWindowSums = []

		#Cycle through all windows of length N, record total overlap in aiWindowSums

		while (fHitEnd == False):
		    if ((iStart+iN) >= len(adAdjWindow)):
				fHitEnd = True
		    dSumCounts = sum(adAdjWindow[iStart:(iStart+iN)])
		    adWindowSums.append(dSumCounts)
		    iStart+=1

		#Find first AminoAcid of best window (lowest total overlap) with length N
		dMin = min(adWindowSums)
		iWinStart = adWindowSums.index(dMin)
		iWinEnd = iWinStart + iN

		bStop = False

		#Add next AA if that does not pus this over the threshhold
		while(bStop==False and iWinEnd< len(dictGenes[key])-1 and len(adAdjWindow[iWinStart:iWinEnd+1]) < iTotLength):
		    if (sum(adAdjWindow[iWinStart:iWinEnd+1])<=iThresh):
				iWinEnd+=1
		    else:
				bStop=True

		#After you've grown the QM as much as possible to the right, try growing to the left.
		#Test the AA before the first AA of the QM, add it if it doesn't put you over iThresh
		while(bStop==False and iWinStart>= 1 and len(aiWindow[iWinStart-1:iWinEnd]) < iTotLength):
		    if (sum(adAdjWindow[iWinStart-1:iWinEnd])<=iThresh):
				iWinStart-=1
		    else:
				bStop=True


		iQuasi = int(sum(adAdjWindow[iWinStart:iWinEnd]))

		#Error Checking

		"""
		if key =="ZP_04318513":
			print "Data"
			print key
			print "Integer Window"
			print aiWindow
			print "Adjusted Window"
			print adAdjWindow
			print "Start,End",iWinStart, iWinEnd
			print adAdjWindow[iWinStart:iWinEnd]
			print dictGenes[key][iWinStart:iWinEnd]
			print "Quasi:",iQuasi
		"""


		tup = (key, dictGenes[key][iWinStart:iWinEnd],iQuasi, iWinStart,iWinEnd,aiWindow[iWinStart:iWinEnd])
		atupQM.append(tup)


		#Error Checking
		"""
		if (key == "VFG1266" or key == "VFG0696" or key=="VFG2059"):
		    print key
		    print dictRefCounts.get(key,"Not in Ref Blast Results")
		    print dictGOICounts.get(key,"Not in Clust Blast Results")
		    print dictGOIGenes[key]
		    print aiWindowSums

		    print iMin, iWinStart

		    print dictGenes[key][iWinStart:iWinStart+iN]
		"""
    return atupQM

######################################################################################
def GetQMOverlap(tupQM,atupHitInfo,fileOut,dictGOIGenes):

		# Convert from Python array index to BLAST position.
		# For example, if the QM comes from the first 20 AA's, it has Start=0 and End=20 for slice notation seq[Start:End]
		# In blast output, it would be Start =1 End = 20

		iMarkerStart = tupQM[3]+1
		iMarkerEnd = tupQM[4]
		iCounter = 1

		# Initialize values used to determine the weights: how much each sequence overlaps the QM (including the original prot)
		iTotalOverlap = 0
		iCurOverlap = 0

		# Print the QM AA's, add it to the overlap total
		fileOut.write(str(tupQM[1]) + '\n')
		#iTotalOverlap += len(tupQM[1])
		iOriginalSeqOverlap = len(tupQM[1])


		iCountOverlappingHits = 0
		for tupHit in atupHitInfo:
			bOverlapsMarker = False

			# tupInfo = strSubId, iQStart, iQEnd, iSubStart,iSubEnd

			strName = tupHit[0]

			# Overlap location on the marker (in Blast format)
			iOLMarkerStart = tupHit[1]
			iOLMarkerEnd = tupHit[2]

	        # Overlap location on the hit (in Blast format)
			iOLHitStart = tupHit[3]
			iOLHitEnd = tupHit[4]

			iShift = iOLMarkerStart - iOLHitStart

			# If the hit is before or after the QM ...
			if (iOLMarkerStart > iMarkerEnd or iOLMarkerEnd < iMarkerStart):
				bOverlapsMarker = False

	        # If the overlapping hit begins within the QM...
			elif (iMarkerStart <= iOLMarkerStart <= iMarkerEnd):
				iOLMarkerEnd = min(iMarkerEnd,iOLMarkerEnd)
				bOverlapsMarker = True

	        # If the overlapping hit ends within the QM...
			elif(iMarkerStart <= iOLMarkerEnd <= iMarkerEnd):
				iOLMarkerStart = max(iMarkerStart,iOLMarkerStart)
				bOverlapsMarker = True

			# If the overlapping hit covers the marker, but doesn't begin or end in the marker, it *completely* overlaps it.
			else:
				iOLMarkerStart = iMarkerStart
				iOLMarkerEnd = iMarkerEnd
				bOverlapsMarker = True

			if (bOverlapsMarker==True):
				if(strName in dictGOIGenes.keys()):

					# The overlapping data from the Marker
					astrOverlapMarker = dictGOIGenes[tupQM[0]][(iMarkerStart):(iMarkerEnd)]

					# The overlapping data from the hit
					astrOverlap = dictGOIGenes[strName][(iOLMarkerStart-iShift-1):(iOLMarkerEnd-iShift)]

					astrOverlap = str("*"*(iOLMarkerStart - iMarkerStart)) + astrOverlap + str("*"*(iMarkerEnd- iOLMarkerEnd))


					# fileOut.write(astrOverlap + ' ' + strName + ' ' + str(iOLMarkerStart) +' '+str(iOLMarkerEnd) + ' '+str(iOLHitStart) + ' ' + str(iOLHitEnd) +' '+  str(tupHit[5]).zfill(3) + '\n')
					fileOut.write(astrOverlap + ' ' + strName + '\n')

				iTotalOverlap+=iOLMarkerEnd-(iOLMarkerStart-1)
				iCountOverlappingHits +=1

		dOriginalSeqWeight = iOriginalSeqOverlap / float(iTotalOverlap + iOriginalSeqOverlap)

		fileOut.write( "Number of Overlapping Seqs: "+ str(iCountOverlappingHits) +'\n')
		fileOut.write( "Sum of Overlapping AA: "+ str(iTotalOverlap) +'\n')
		fileOut.write( "Original Sequence's Weight: "+ "{:3.2f}".format(dOriginalSeqWeight) +'\n' + '\n')

		return iTotalOverlap

####################################################################################


def UpdateQMHeader(atupQM,dictGOIHits,dictRefHits,strQMOut,dictGOIGenes,bUpdateHeader=False):
    #Each QM_tuple has the values: (name, window, overlapvalue,iMarkerStart,iMarkerEnd,aOverlap)


	atupUpdatedQM = []

	with open(strQMOut, 'w') as fQMOut:
		for tupQM in atupQM:
			fQMOut.write("***********************************************************************" + '\n')

			# Write the QM out
			fQMOut.write(">" + str(tupQM[0]) + "_QM" + str(tupQM[2]) + "_#" + '\n')
			fQMOut.write(str(tupQM[1]) + '\n' + '\n')

			fQMOut.write("Overlap array for QM (before ^(1/4) transformation):" + '\n' + str(tupQM[5]) + '\n' + '\n')

			fQMOut.write("Overlap among GOI Seqs:" + '\n')
			iGOIOverlap = GetQMOverlap(tupQM,dictGOIHits[tupQM[0]],fQMOut,dictGOIGenes)


			if tupQM[0] in dictRefHits.keys():
				fQMOut.write("Overlap among Ref Seqs:" + '\n')
				iRefOverlap = GetQMOverlap(tupQM,dictRefHits[tupQM[0]],fQMOut,dictGOIGenes)
			else:
				iRefOverlap = 0


			dFinalWeight = len(tupQM[1]) / float(iGOIOverlap + iRefOverlap + len(tupQM[1]))
			strFinalWeight = "{:3.2f}".format(dFinalWeight)
			fQMOut.write("Final weight: " + strFinalWeight + '\n\n')
			aQM = list(tupQM) + [strFinalWeight]
			atupUpdatedQM.append(aQM)

	return atupUpdatedQM

###############################################################################
def PrintQuasiMarkers(atupQM, fileOut):
	iCounter = 0
 	strName = ""


	for tup in atupQM:
		if str(tup[0]) != strName:
			iCounter =1
		else:
			iCounter+=1
		fileOut.write(">" + str(tup[0]) + "_QM" + str(tup[2]) + "_#" +str(iCounter).zfill(2) + "_w"  + '\n')
		fileOut.write(str(tup[1]) + '\n')
		#fileOut.write(str(tup))
		strName = str(tup[0])


	return
