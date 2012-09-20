#!/usr/bin/python


#  cd /media/OS/thesis/src/
#python process_blast.py < /media/OS/thesis/tmp/clustTAB30.blast > /media/OS/thesis/tmp/out.txt

import re
import sys
import csv
import argparse

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Remove hits with indetity >=i and length <=l from sequence.')
parser.add_argument('--i',default = .90, type=float, dest='dID', help='Enter the identity cutoff. Examples: .90, .85, .10,...')
parser.add_argument('--l', default = .10, type=float, dest='dL', help='Enter maximum length for region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
parser.add_argument('--tabs', default = False, type=bool, dest='fTabs', help='Set to True if you would like to print: GeneName, #AAs after process, ##AA initially ')
parser.add_argument('--fasta', type=file, dest='fFile', help='Enter the path and name of your fasta file.')
parser.add_argument('--abs_l', default=0, type=int, dest='iRegionLength', help='Enter an integer here to remove regions of X length or greater. This will cause the program to ignore any other parameters for length')
args = parser.parse_args()
             
"""
Add a log file
-Number of genes in input file
-Number of matches in blast file
-Number of unique query genes in blast file
-Number of genes that did not drop any domains
-Number of genes that lost more than 90% of data

"""             
             
dIDcutoff = args.dID * 100
dLengthcutoff = args.dL

###########################################################
def getGeneData ( fileFasta):
    dictGeneData = {}
    
    for gene in SeqIO.parse(fileFasta, "fasta"):
            dictGeneData.setdefault(gene.id.strip(), str(gene.seq))
        
    
    if dictGeneData.has_key(''):        
        del dictGeneData[''] 
    return dictGeneData

##############################################################################
#Have a dictionary that is (gene, length)

dictGeneData = getGeneData(args.fFile)


##############################################################################
#Read in the blast output line by line
#When the program finds a new QueryGene:
    #Add the last gene's abWindow to dictGeneWindows -- (gene, abWindow) 
    #Make a new abWindow
#When the program finds the same QueryGene:
    #Make a new abWindow the length of the gene +1, set all of it equal to True
    #If the region in the blast hit satisfies our length and ID parameters
        #Set the corresponding location in abWindow equal to FALSE

#Add the last gene when you finish the file
#Delete the 


strCurQuery = ""
dictGeneWindows = {}
abWindow =[]
iGeneCount = 0

for aLine in csv.reader( sys.stdin, csv.excel_tab ):
    strQueryID = aLine[0]
    strSubId = aLine[1]
    dIdentity =float(aLine[2])
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
    
    if dictGeneData.has_key(strQueryID.strip()):
        if strQueryID != strCurQuery:
            dictGeneWindows.setdefault(strCurQuery, abWindow)        
            iGeneCount = iGeneCount+1        
            
            strCurQuery = strQueryID
                  
            
            
            # we will ignore abWindow[0]
            abWindow = []
            for i in range(iQLength+1):
                abWindow.append(True)
        
        
        dMatchLength = (iAln) / float(iQLength)
        abWindow[0] = False

        #If user gave "abs_l" parameter, use that to determine what regions to eliminate
        if (args.iRegionLength>0):
            if (dIdentity >= dIDcutoff) and (iAln >= args.iRegionLength) and (strQueryID!=strSubId):
                for i in range(iQStart, iQEnd+1):
                    abWindow[i]=False
            
        #Else: Mask high-identity, low-length regions using (alignment length / query length)
        elif (dIdentity >= dIDcutoff) and (dMatchLength <= dLengthcutoff) and (strQueryID!=strSubId):
            for i in range(iQStart, iQEnd+1):
                abWindow[i]=False

#Add the last gene
dictGeneWindows.setdefault(strCurQuery.strip(), abWindow)        
iGeneCount = iGeneCount+1 

#############################################################################
#Print iGeneCount if fTabs is True. This is for debugging
del dictGeneWindows[""]
if args.fTabs == True:
    for key in dictGeneWindows:
                print key,",", dictGeneWindows[key].count(True),",", len(dictGeneWindows[key]), ",",args.dID,",", args.dL



###########################################################################
#Replace domains with X's
for key in dictGeneData:
        if dictGeneWindows.has_key(key):
            abWindow = dictGeneWindows[key]
            strGene = list(dictGeneData[key])
            #Blast results start at 1, and Python starts at 0. Deleting abWindow[0] lines them up
            del abWindow[0]
            for i in range(0,len(abWindow)):
                if abWindow[i] == False:
                    strGene[i] = "X"
            strGene = "".join(strGene)
            dictGeneData[key] =strGene


##############################################################################
#Print out the genes
strGeneName = ""
iCount = 0

     
for key in dictGeneData:
    if dictGeneWindows.has_key(key):
        strGeneName = ">" + key
    else: 
        strGeneName = ">" + key + " NO BLAST MATCHES"
        
    
    if args.fTabs == False:
        print strGeneName
        print re.sub("(.{80})","\\1\n",dictGeneData[key],re.DOTALL)
    iCount = iCount+1
