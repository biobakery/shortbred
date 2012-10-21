#!/usr/bin/env python


#This is a modified version of process_blast.py. That program used "bool" values
#to determine if an amino acid was X'ed out. This program instead records the count
#of regions that overlapped a given AA.

#****************************************************************************
#--Example Use--
#  python process_blast_int.py --fasta /home/jim/tmp/clust.faa --goi /home/jim/tmp/blast_clust_to_clust.txt --ref /home/jim/tmp/blast_clust_to_clust.txt


import re
import sys
import csv
import argparse
import os
import subprocess
import glob

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

"""
parser = argparse.ArgumentParser(description='Remove hits with indetity >=i and length <=l from sequence.')
parser.add_argument('--id',default = .90, type=float, dest='dID', help='Enter the identity cutoff. Examples: .90, .85, .10,...')
parser.add_argument('--len', default = .10, type=float, dest='dL', help='Enter maximum length for region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
parser.add_argument('--tabs', default = False, type=bool, dest='fTabs', help='Set to True if you would like to print: GeneName, #AAs after process, ##AA initially ')
parser.add_argument('--fasta', type=file, dest='fGOIFasta', help='Enter the path and name of your fasta file.')
parser.add_argument('--abs_l', default=0, type=int, dest='iRegionLength', help='Enter an integer here to remove regions of X length or greater. This will cause the program to ignore any other parameters for length')
parser.add_argument('--ref', type=file, dest='fRefBlast', help='Enter the path and name of the blast results from the refrence db.')
parser.add_argument('--goi', type=file, dest='fGOIBlast', help='Enter the path and name of the blast results from the goi db.')
parser.add_argument('--ko', type=str, dest='sKO', help='Enter \"ref\",\"goi\" or \"both\".')
parser.add_argument('--winlength', type=int, dest='iWinLength', help='Enter window length')

args = parser.parse_args()
"""             
"""
Add a log file
-Number of genes in input file
-Number of matches in blast file
-Number of unique query genes in blast file
-Number of genes that did not drop any domains
-Number of genes that lost more than 90% of data
             


dIDcutoff = args.dID * 100
dLengthcutoff = args.dL

"""
###############################################################################
def MakeFamilyFastaFiles ( dictFams, fileFasta, dirOut):
    
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
def ClusterFams(dirClust):
    
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
        fileUC    = dirUC + os.sep + os.path.basename(fileFasta) + ".uc"
        subprocess.check_call(["usearch6", "--cluster_fast", str(fileFasta), "--uc", str(fileUC), "--id", ".95","--centroids", str(fileClust)])

    
    ageneAllGenes = []

    for fileFasta in glob.glob(dirCentroids+os.sep+'*.faa'):
        for gene in SeqIO.parse(fileFasta, "fasta"):
            ageneAllGenes.append(gene)
            
    
    SeqIO.write(ageneAllGenes, dirClust + os.sep + "clust.faa", "fasta")

###############################################################################
def printMap(strUCMap,strTxtMap):
   
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
    dictGeneData = {}
    
    for gene in SeqIO.parse(fileFasta, "fasta"):
            dictGeneData.setdefault(gene.id.strip(), str(gene.seq))
        
    
    if dictGeneData.has_key(''):        
        del dictGeneData[''] 
    return dictGeneData

##############################################################################
#Make a dictionary of form (Gene, [0,0,0,0,1,1...]), where the number indicates
#the number of times an amino acid overlaps with a region in the blast output

#Read in the blast output line by line
#When the program finds a new QueryGene:
    #Add the last gene's aiCounts to dictGeneWindows -- (gene, abWindow) 
    #Make a new aiCounts
#When the program finds the same QueryGene:
    #If the region in the blast hit satisfies our length and ID parameters
        #Set the corresponding location in abWindow equal to (count+1)

#Add the last gene when you finish the file

def getOverlapCounts (fileBlast, dIDcutoff, dLengthMin, dLengthcutoff, iOffset, iRegionLength, dictFams={}):

    strCurQuery = ""
    dictAAOverlapCounts = {}
    aiCounts =[]
    iGeneCount = 0
    
    for aLine in csv.reader( open(fileBlast), csv.excel_tab ):
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

        if strQueryID != strCurQuery:
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

        if (len(dictFams) > 0):
            bNotSameFam = (dictFams.get(strQueryID,"NoQuery")!=dictFams.get(strSubId,"NoSubject"))
        #Mask high-identity, low-length regions using (alignment length / query length)
        else:
            bNotSameFam = True
        
        if (dIdentity >= dIDcutoff) and (dMatchLength <= dLengthcutoff) and (strQueryID!=strSubId) and (dMatchLength >= dLengthMin) and bNotSameFam:
            for i in range(iQStart-1+iOffset, iQEnd-iOffset):
                aiCounts[i]=aiCounts[i]+1
    
    #Once the loop is done, remember to add the last window.
    dictAAOverlapCounts.setdefault(strCurQuery.strip(), aiCounts)        
    iGeneCount = iGeneCount+1 
    
    del dictAAOverlapCounts[""]

    return dictAAOverlapCounts

###########################################################################
#Take the genes in setNames, look in dictKnockOut to see if they have a 
#region of length N without any overlap. Returns set of genes with markers.


def CheckForMarkers(setGenes, dictKnockOut, iN):

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
def CheckForQuasiMarkers(setGenes, dictKnockOut, dictGenes, iN):
    
    #Only run this on the leftover genes
    # "n" = minimum window length
    #For each one, sum up the values from [0:n], then [1:n+1]...
    #Store these in an array of length (len(gene)-n)
    #FInd the minimum value in this array
    #Take its index
    #Your window is [index:index+n]
    
    #Get the appropriate string from dictGOIgenes    
    #add it to dictQM
    
    
    #Return dictQM with these windows

    #tuple with 3 values (name, window, overlapvalue)
    
    
    atupQM = []

    
    
    for key in setGenes:
        aiWindow = dictKnockOut[key]
        iStart = 0
        iMin = 0
        iSumCounts = 0
        fHitEnd = False
        aiWindowSums = []
 
        

        
        #Cycle through all windows of length N, record total overlap
        
        while (fHitEnd == False):
            if ((iStart+iN) >= len(aiWindow)):
                fHitEnd = True                
            iSumCounts = sum(aiWindow[iStart:(iStart+iN)])        
            aiWindowSums.append(iSumCounts)
            iStart+=1
            
        #Find first AminoAcid of best window (lowest total overlap) with length N
        iMin = min(aiWindowSums)
        iWinStart = aiWindowSums.index(iMin)            
        iWinEnd = iWinStart + iN
        
        bStop = False        
        
        #If the next AA has (overlap==0), then extend the window to include it                        
        while(bStop==False and iWinEnd< len(dictGenes[key])-1):
            if (dictGenes[key][iWinEnd+1]==0):
                iWinEnd+=1
                print key, "this ran"
            else:
                bStop=True
        
        tup = (key, dictGenes[key][iWinStart:iWinEnd],iMin)       
        atupQM.append(tup)
        
        #In the knockouty dictionary, set the QM region you just took to have "9999"
        #for each AA. Then it wont be used again for the second QM window.
        dictKnockOut[key][iWinStart:iWinEnd] = [9999]*(iWinEnd-iWinStart +1)
        
        
        
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

"""
##############################################################################
#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
dictGOIGenes = getGeneData(args.fGOIFasta)
dictRefCounts = getOverlapCounts(args.fRefBlast)
dictGOICounts = getOverlapCounts(args.fGOIBlast)



#If a gene has 0 valid hits in the ref database, make an array of 0's
#so the program knows that nothing overlapped with the gene
setGOINotInRef = set(dictGOIGenes.keys()).difference(set(dictRefCounts.keys()))

if len(setGOINotInRef)>0:
    for sGene in setGOINotInRef:
        dictRefCounts[sGene] = [0]*len(dictGOIGenes[sGene])

#If a gene has 0 valid hits in the GOI database (unlikely), make an 
#array of 0's so the program knows that nothing overlapped with the gene    

setGOINoHits = set(dictGOIGenes.keys()).difference(set(dictGOICounts.keys()))

if len(setGOINoHits)>0:
    for sGene in setGOINoHits:
        dictGOICounts[sGene] = [0]*len(dictGOIGenes[sGene])


#Get dict of counts for (Ref+GOI)
dictBoth = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys())) 

for sGene in setRefGOI:
    aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]))]
    dictBoth[sGene] = aiSum
    


###########################################################################
#Check for genes that have "marker windows": windows of length N that do not
#overlap with anything in the "overlap reference"


setHasMarkers = CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictBoth.keys()), dictBoth, args.iWinLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)
"""
#Removed this code, was used to make class markers

"""
setHasClassMarkers = CheckForMarkers(setLeftover.intersection(dictRefCounts.keys()), dictRefCounts, args.iWinLength)
print "Genes with Class Markers:",  len(setHasClassMarkers)


setLeftover = setLeftover.difference(setHasClassMarkers)
"""
"""
dictQuasiMarkers = CheckForQuasiMarkers(setLeftover, dictBoth, dictGOIGenes,args.iWinLength)



###########################################################################
#Replace AA's with X's in True Markers
for key in setHasMarkers:
        if dictBoth.has_key(key):
            aiWindow = dictBoth[key]
            strGene = list(dictGOIGenes[key])

            for i in range(0,len(aiWindow)):
                if aiWindow[i] >= 1:
                    strGene[i] = "X"
            strGene = "".join(strGene)
            dictGOIGenes[key] =strGene

###########################################################################
#Replace AA's with X's in Class Markers
for key in setHasClassMarkers:
        if dictRefCounts.has_key(key):
            aiWindow = dictRefCounts[key]
            strGene = list(dictGOIGenes[key])

            for i in range(0,len(aiWindow)):
                if aiWindow[i] >= 1:
                    strGene[i] = "X"
            strGene = "".join(strGene)
            dictGOIGenes[key] =strGene
           
###########################################################################
#Add in the QuasiMarkers

for key in dictQuasiMarkers:
    dictGOIGenes[key]=dictQuasiMarkers[key][0]
   
##############################################################################
#Print out the genes
strGeneName = ""
iCount = 0


#f.write( "Database, ClustID, PctID, PctLength, WinLength, Quasi-Threshhold Used, Num Genes After Clust, # TM's, # QM's, Windows, Num Genes Missing Windows, Sum Overlap for QM's")
print("VF,95%," + str(args.dID) + "," + str(args.dL) + "," + str(args.iWinLength) + ",999," + str(len(dictGOIGenes)))  


print "GOI Genes:", len(dictGOIGenes)
print "Ref Valid Hits:", len(dictRefCounts)
print "GOI Valid Hits:",len(dictGOICounts)

print ""
print "Window Length: ", args.iWinLength

print "Genes with True Markers:", len(setHasMarkers)
print "Genes with Quasi-Markers:", len(dictQuasiMarkers)
print "Union of all Genes with Markers:", len(setHasMarkers.union(dictQuasiMarkers.keys()))

print "Length of dictGOIGenes:",len(dictGOIGenes)
     
for key in dictGOIGenes:
    if key in setHasMarkers:
        strGeneName = ">" + key + "_TM"
    #elif key in setHasClassMarkers:
     #   strGeneName = ">" + key + "_CM"
    elif key in dictQuasiMarkers:
        strGeneName = ">" + key + "_QM" + str(dictQuasiMarkers[key][1])
    else:
        strGeneName = ">" + key + "_OTH"
         
    
    print strGeneName
    print re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)
    iCount = iCount+1

"""