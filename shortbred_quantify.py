#!/usr/bin/env python

import sys
import argparse
import subprocess
import csv
import re
import os
import datetime
import shutil

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers and wgs file as input, and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--wgs', type=str, dest='sWGS', help='Enter the path and name of the genes of interest file (proteins).')

#Output
parser.add_argument('--SBhits', type=str, dest='sHits', help='ShortBRED will print the hits it considers positives to this file.', default="")
parser.add_argument('--blastout', type=str, dest='strBlast', default="out.blast",help='Enter the path and name of the blastoutput.')

#Parameters
parser.add_argument('--id', type=float, dest='dID', help='Enter the percent identity for the match', default = .95)
parser.add_argument('--tmid', type=float, dest='dTMID', help='Enter the percent identity for a TM match', default = .95)
parser.add_argument('--qmid', type=float, dest='dQMID', help='Enter the percent identity for a QM match', default = .90)
parser.add_argument('--alnlength', type=int, dest='iAlnLength', help='Enter the minimum alignment length. The default is 20', default = 20)
parser.add_argument('--alnTM', type=int, dest='iAlnMax', help='Enter a bound for TM alignments, such that aln must be>= min(markerlength,alnTM)', default = 25)


parser.add_argument('--tmp', type=str, dest='sTmp', default =os.getcwd() +os.sep + "tmp",help='Enter the path and name of the tmp directory.')
parser.add_argument('--length', type=int, dest='iLength', help='Enter the minimum length of the markers.')
parser.add_argument('--threads', type=int, dest='iThreads', help='Enter the number of CPUs available for usearch.', default=1)
parser.add_argument('--notmarkers', type=str, dest='strNM',default="N", help='.')

#DB Note - Maybe ask Nicola how to remove usearch6 output

 

args = parser.parse_args()

dirTmp = args.sTmp
if not os.path.exists(dirTmp):
    os.makedirs(dirTmp)

    


strBase = os.path.splitext(os.path.basename(args.sMarkers))[0]

if(args.sHits==""):
    strHitsFile = args.sTmp + os.sep + "SBhits.txt"
else:
    strHitsFile = args.sHits

log = open(str(dirTmp + os.sep + strBase+ ".log"), "w")
log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n SEARCH PARAMETERS \n")
log.write("Match ID:" + str(args.dID) + "\n")
log.write("Alignment Length:" + str(args.iAlnLength) + "\n")
log.write("TM id:" + str(args.dTMID) + "\n")
log.write("QM id:" + str(args.dQMID) + "\n")
log.close()

###############################################################################
# SUM TOTAL MARKER LENGTH FOR EACH PROT FAMILY

dictMarkerLen = {}
dictMarkerLenAll = {}


#Load the WGS data. For each gene stub, copy in all the matching genes
for seq in SeqIO.parse(args.sMarkers, "fasta"):
	if args.strNM=="N":
		mtchStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',seq.id)
		strStub = mtchStub.group(1)
	else:
		strStub = seq.id
	dictMarkerLenAll[strStub] = len(seq) + dictMarkerLenAll.get(strStub,0)   
 	dictMarkerLen[seq.id] = len(seq)
    

###############################################################################
#USE USEARCH, CHECK WGS NUCS AGAINST MARKER DB

#Make a database from the markers

strDBName = str(args.sMarkers) + ".udb"
strSearchResults = args.sMarkers +".blast"

#Note: Cannot limit threads used in creating of usearch database
p = subprocess.check_call(["usearch6", "--makeudb_usearch", str(args.sMarkers), "--output", strDBName])


if args.strNM=="N":
    #Use usearch to check for hits (usearch local)
    subprocess.check_call(["time","-o", args.sMarkers+ ".time","usearch6", "--usearch_local", str(args.sWGS), "--db", str(strDBName), "--id", str(args.dID),"--blast6out", args.strBlast,"--threads", str(args.iThreads)])
else:
    #Use usearch to check for hits (usearch local)
    subprocess.check_call(["time","-o",args.sMarkers +".time","usearch6", "--usearch_local", str(args.sWGS), "--db", str(strDBName), "--id", ".90","--blast6out", args.strBlast,"--threads", str(args.iThreads)])



fileHits = open(strHitsFile,'w') 

#Go through the blast hits, for each prot family, print out the number of hits
dictBLAST = {}    
for aLine in csv.reader( open(args.strBlast), csv.excel_tab ):
    mtchTM = re.search(r'_TM.*',aLine[1])
    if (mtchTM):
        dID = args.dTMID
        iAln = min(dictMarkerLen[aLine[1]] ,args.iAlnMax)
    else:
        dID = args.dQMID
        iAln = args.iAlnLength
    
    if args.strNM=="N":
        mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])    
        strProtFamily = mtchProtStub.group(1)
        if (int(aLine[3])>= iAln and (float(aLine[2])/100.0) >= dID):
            dictBLAST.setdefault(strProtFamily,set()).add((aLine[0]))
            fileHits.write('\t'.join(aLine) + '\n')
    else:
        strProtFamily = aLine[1]
	if (int(aLine[3])>= iAln):
            dictBLAST.setdefault(strProtFamily,set()).add((aLine[0]))
            fileHits.write('\t'.join(aLine) + '\n')        		

fileHits.close()                
	
    
for strProt in dictBLAST.keys():
    print strProt + "\t" +  str(float(len(dictBLAST[strProt]))/dictMarkerLenAll[strProt])
    
