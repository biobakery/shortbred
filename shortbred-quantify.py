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


parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name of the genes of interest file (proteins).')

parser.add_argument('--wgs', type=str, dest='sWGS', help='Enter the path and name of the genes of interest file (proteins).')

parser.add_argument('--tmp', type=str, dest='sTmp', help='Enter the path and name of the tmp directory.')
parser.add_argument('--length', type=int, dest='iLength', help='Enter the minimum length of the markers.')


args = parser.parse_args()

###############################################################################
# SUM TOTAL MARKER LENGTH FOR EACH PROT FAMILY

dictMarkerLen = {}

#Load the WGS data. For each gene stub, copy in all the matching genes
for seq in SeqIO.parse(args.sMarkers, "fasta"):
    mtchStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',seq.id)
    dictMarkerLen[mtchStub.group(1)] = len(seq) + dictMarkerLen.get(mtchStub.group(1),0)   
    

###############################################################################
#USE USEARCH, CHECK WGS NUCS AGAINST MARKER DB

#Make a database from the markers

strDBName = args.sMarkers + ".udb"
strSearchResults = args.sMarkers +".blast"

p = subprocess.Popen(["usearch6", "--makeudb_usearch", args.sMarkers, "--output", strDBName],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
out, err = p.communicate()

dQueryCov = (float(args.iLength)*3)/float(100)
print dQueryCov
dQueryCov = round(dQueryCov,2)
strQC = str(dQueryCov)
print "The value is: ", strQC

#Use usearch to check for hits (usearch global)
p = subprocess.Popen(["usearch6", "--usearch_local", args.sWGS, "--db", strDBName, "--id", "0.90", "--query_cov", strQC, "--blast6out", strSearchResults])





#Use usearch to checkj for hits (usearch local)

iReads = 0
iHits = 0

for line in open(args.sWGS):
    mtch = re.search(r'^\>',line)
    if (mtch):
        iReads+=1

#Go through the blast hits, for each prot family, print out the number of hits
dictBLAST = {}    
for aLine in csv.reader( open( strSearchResults), csv.excel_tab ):
    #print aLine[1]
    mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])    
    strProtFamily = mtchProtStub.group(1)        
    dictBLAST.setdefault(strProtFamily,set()).add((aLine[0]))
    if (mtchProtStub):
        iHits+=1
            
    
for strProt in dictBLAST.keys():
    print strProt + "\t" +  str(float(len(dictBLAST[strProt]))/dictMarkerLen[strProt])
    

print "Reads in WGS file:" + " " + str(iReads)
print "Total Reads Matched:" + " " + str(iHits)
print "Share of WGS matched:" + " " + str(round(float(iHits)/iReads,5))