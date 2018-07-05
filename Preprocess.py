#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocess inputfile to be used by ShortBRED
===============================================
* Loads proteins sequences from FASTA
* Check whether the input file can be used by ShortBRED directly
* If it cannot be used, create a modified version 
* Replace all illegal characters in the protein ids with underscores
* Replace the strings '_TM','_QM',and '_JM' to the strings '_tm','_tm',and '_jm'
* Make sure the protein ids are unique
===============================================
Author: Jingjing Tang (jatangne@gmail.com)
"""

import re
import os
import src
import sys
import argparse

try:
    import Bio
except ImportError:
    print("\nShortBRED was unable to load Biopython. Please check to make sure you have Biopython installed (http://biopython.org/wiki/Main_Page), and that its directory is in your PYTHONPATH. \n")
    sys.exit(1)

from Bio import SeqIO


def CheckFastaForBadProtNames(fileFasta):
    reBadChars=re.compile(r'[\\\/\*\=\:\'\[\]\.\;\,]')
    reMarkerTypes = re.compile(r'_([TJQ]M)')
    setProtNames = set()
    
    for gene in SeqIO.parse(fileFasta, "fasta"):
        mtchBad = reBadChars.search(gene.id)
        mtchMarkerType = reMarkerTypes.search(gene.id)
        
        if mtchBad != None or mtchMarkerType != None or gene.id in setProtNames:
            return False
        else:
            setProtNames.add(gene.id)
    return True


def preprocessInput(inputfile, outputfile):    
    fasta_sequences = SeqIO.parse(open(inputfile),'fasta')
    dictProtNames = {}
    with open(outputfile, 'w') as ofile:
        for gene in fasta_sequences:
            gene.id = re.sub('[^0-9a-zA-Z]+', '_', gene.id)
            gene.id = re.sub('_TM', '_tm',gene.id)
            gene.id = re.sub('_QM', '_qm',gene.id)
            gene.id = re.sub('_JM', '_jm',gene.id)
            if gene.id in dictProtNames.keys():
                dictProtNames[gene.id] += 1
                gene.id = gene.id + "_" + str(dictProtNames[gene.id])
            else:
                dictProtNames[gene.id] = 1
            SeqIO.write(gene, ofile, 'fasta')
    return

##############################################################################

parser = argparse.ArgumentParser( formatter_class=argparse.ArgumentDefaultsHelpFormatter )
parser.add_argument( "--input", type = str, dest="inputfile", help="Enter the path and name of the proteins of interest file." )
args = parser.parse_args()

checkresult = CheckFastaForBadProtNames(args.inputfile)

if checkresult:
    sys.stderr.write( "Processing complete! The proteins of interest file can be used by ShortBRED directly\n")
else:
    pathlist = args.inputfile.split("/")
    filename = "Modified_" + pathlist[-1]
    dirTmp = "/".join(pathlist[:-1])
    outputfile = dirTmp + '/' + filename
#    if os.path.isfile(outputfile):
#        os.remove(outputfile)
    preprocessInput(args.inputfile, outputfile)
    sys.stderr.write( "Processing complete! The proteins of interest file cannot be used by ShortBRED directly\n")
    sys.stderr.write( "\nProcessing complete! New FASTA file saved to " + outputfile + "\n")
    sys.stderr.write( "\nNOTE: Please use the modified proteins of interest file instead.\n\n")

