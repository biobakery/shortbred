#!/usr/bin/env python

import sys
import argparse
import subprocess
import csv
import re
import os

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')

#Input
parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--wgs', type=str, dest='sWGS', help='Enter the path and name of the genes of interest file (proteins).')

#Output
parser.add_argument('--results', type=str, dest='sResults', help='Enter the path and name of the results file.')

#Parameters
parser.add_argument('--tmp', type=str, dest='sTmp', help='Enter the path and name of the tmp directory.')
parser.add_argument('--length', type=int, dest='iLength', help='Enter the minimum length of the markers.')
parser.add_argument('--threads', type=int, dest='iThreads', help='Enter the number of CPUs available for usearch.')
parser.add_argument('--notmarkers', type=str, dest='strNM',default="N", help='.')


args = parser.parse_args()

###############################################################################
# SUM TOTAL MARKER LENGTH FOR EACH PROT FAMILY

dictMarkerLen = {}

#Load the WGS data. For each gene stub, copy in all the matching genes
for seq in SeqIO.parse(args.sMarkers, "fasta"):
	if args.strNM=="N":
		mtchStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',seq.id)
		strStub = mtchStub.group(1)
	else:
		strStub = seq.id
	dictMarkerLen[strStub] = len(seq) + dictMarkerLen.get(strStub,0)   
    

###############################################################################
#USE USEARCH, CHECK WGS NUCS AGAINST MARKER DB

#Make a database from the markers

strDBName = args.sMarkers + ".udb"
strSearchResults = args.sMarkers +".blast"

#Note: Cannot limit threads used in creating of usearch database
p = subprocess.check_call(["usearch6", "--makeudb_usearch", args.sMarkers, "--output", strDBName])

dQueryCov = (float(args.iLength)*3)/float(100)
print dQueryCov
dQueryCov = round(dQueryCov,2)
strQC = str(dQueryCov)
print "The value is: ", strQC

#Use usearch to check for hits (usearch local)
subprocess.check_call(["usearch6", "--usearch_local", args.sWGS, "--db", strDBName, "--id", "0.95", "--blast6out", strSearchResults,"--threads", str(args.iThreads), "--gapopen", "20"])


#Use usearch to checkj for hits (usearch local)



#Go through the blast hits, for each prot family, print out the number of hits
dictBLAST = {}    
for aLine in csv.reader( open( strSearchResults), csv.excel_tab ):
	#print aLine[1]
	if args.strNM=="N":
		mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])    
		strProtFamily = mtchProtStub.group(1)
	else:
		strProtFamily = aLine[1]
	dictBLAST.setdefault(strProtFamily,set()).add((aLine[0]))
	


    
for strProt in dictBLAST.keys():
    print strProt + "\t" +  str(float(len(dictBLAST[strProt]))/dictMarkerLen[strProt])
    
