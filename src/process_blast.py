#!/usr/bin/env python

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
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

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
def ClusterFams(dirClust, dCLustID, strOutputFile):
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

        iSeqCount = 0        
        #Count seqs, if more than one, then align them        
        for seq in SeqIO.parse(fileFasta, "fasta"):
            iSeqCount+=1
            
        
        if iSeqCount>1:
            #Call muscle to produce an alignment
            subprocess.check_call(["muscle", "-in", str(fileFasta), "-out", str(fileAlign)])          
            
            #Call cons or em_cons from the EMBOSS package to produce a consensus sequence
            subprocess.check_call(["cons", "-seq", str(fileAlign), "-outseq", str(fileClust)])
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
def getOverlapCounts (fileBlast, dIDcutoff, dLengthMin, dLengthcutoff, iOffset, iRegionLength, dictFams={}):
#Makes a dictionary of form (GeneID, [0,0,0,0,1,1...]), where the number indicates
#the number of times an amino acid overlaps with a region in the blast output.

#Read in the blast output line by line
#When the program finds a new QueryGene:
    #Add the last gene's aiCounts to dictGeneWindows -- (gene, abWindow) 
    #Make a new aiCounts
#When the program finds the same QueryGene:
    #If the region in the blast hit satisfies our length and ID parameters
        #Set the corresponding location in abWindow equal to (count+1)

#Add the last gene when you finish the file

    strCurQuery = ""
    dictAAOverlapCounts = {}
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
        #DB Note - Edited so we do not add the blank strCurQuery on the first line.
        if strQueryID != strCurQuery:
            if iLine>1:
                dictAAOverlapCounts.setdefault(strCurQuery, aiCounts)        
            iGeneCount = iGeneCount+1        
            strCurQuery = strQueryID
                  
            
            aiCounts = []
            for i in range(iQLength):
                aiCounts.append(0)
        
        dMatchLength = (iAln) / float(iQLength)

        """
        #If user gave "abs_l" parameter, use that to determine what regions to eliminate
        if (iRegionLength>0):
            if (dIdentity >= dIDcutoff) and (iAln >= iRegionLength) and (strQueryID!=strSubId):
                #(Blast starts count at 1, but our array starts at 0, so we subtract 1. )                
                for i in range(iQStart-1, iQEnd):
                    aiCounts[i]=aiCounts[i]+1
        """
        
        #DB Note - This code was part of an old idea, to not knockout seqs that belong to the same family if the user pre-specified families.  
        """
        if (len(dictFams) > 0):
            bNotSameFam = (dictFams.get(strQueryID,"NoQuery")!=dictFams.get(strSubId,"NoSubject"))
        #Mask high-identity, low-length regions using (alignment length / query length)
        else:
            bNotSameFam = True
        """
        bNotSameFam = True
        
        if (dIdentity >= dIDcutoff) and (dMatchLength <= dLengthcutoff) and (strQueryID!=strSubId) and (dMatchLength >= dLengthMin) and bNotSameFam:
            for i in range(iQStart-1+iOffset, iQEnd-iOffset):
                aiCounts[i]=aiCounts[i]+1
    
    #Once the loop is done, remember to add the last window.
    dictAAOverlapCounts.setdefault(strCurQuery.strip(), aiCounts)        
    iGeneCount = iGeneCount+1 

    return dictAAOverlapCounts

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

    #A QM_tuple has 3 values (name, markers, overlapvalue)
    
    
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
        
        tup = (key, dictGenes[key][iWinStart:iWinEnd],iQuasi)       
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

###############################################################################
def PrintQuasiMarkers(atupQM, fileOut):
    iCounter = 0
    strName = ""
    
    
    for tup in atupQM:
        if str(tup[0]) != strName:
            iCounter =1  
        else:
            iCounter+=1
        fileOut.write(">" + str(tup[0]) + "_QM" + str(tup[2]) + "_#" +str(iCounter).zfill(2) + '\n')
        fileOut.write(str(tup[1]) + '\n')
        strName = str(tup[0])
    
        
    return
