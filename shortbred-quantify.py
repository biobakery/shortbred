#!/usr/bin/env python

import sys
import argparse
import subprocess
import csv
import re

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name of the genes of interest file (proteins).')

parser.add_argument('--wgs', type=str, dest='sWGS', help='Enter the path and name of the genes of interest file (proteins).')

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
strSearchResults = "results.blast"

subprocess.check_call(["usearch6", "--makeudb_usearch", args.sMarkers, "--output", strDBName])



#Use usearch to checkj for hits (usearch global)
subprocess.check_call(["usearch6", "--usearch_global", args.sWGS, "--db", strDBName, "--id", "0.9", "--blast6out", strSearchResults])




#Use usearch to checkj for hits (usearch local)




#Go through the blast hits, for each prot family, print out the number of hits
dictBLAST = {}    
for aLine in csv.reader( open(strSearchResults), csv.excel_tab ):
    mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',aLine[1])    
    strProtFamily = mtchProtStub.group(1)        
    dictBLAST.setdefault(strProtFamily,set()).add((aLine[0]))
            
    
for strProt in dictBLAST.keys():
    print strProt, float(len(dictBLAST[strProt]))/dictMarkerLen[strProt]
    


