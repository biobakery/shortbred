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


import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO


from collections import Counter




import re
import sys
import argparse


parser = argparse.ArgumentParser(description='Compares file of windows against orginial, clustered file to count windows, check for missing windows.')
parser.add_argument('--cf', type=file, dest='fClusteredFile', help='Enter the clustered file from the project.')
parser.add_argument('--list', type=bool,default=False, dest='bList', help='Set to \'True\' if you want the gene names.')
args = parser.parse_args()



def getSeqs(fileFasta):

    agSeqs = []

    for seq in SeqIO.parse(fileFasta, "fasta"):
        agSeqs.append(seq)

    return agSeqs

agClustSeqs = getSeqs(args.fClusteredFile)
agWindowSeqs = getSeqs(sys.stdin)

##############################################################################
#Prep Data

setHasTM = set()
setHasQM = set()
iTM = 0
iQM = 0


astrWindowNames = []
for i in agWindowSeqs:
    mtchName = re.search(r'(.*)(\_[TCQ]M[0-9]*)(\_\#[0-9]*)',i.id)
    strName = mtchName.group(1)
    strMarkerType = mtchName.group(2)
    if strMarkerType == "_TM":
        setHasTM.add(strName)
        iTM+=1
    if strMarkerType[:3] == "_QM":
        setHasQM.add(strName)
        iQM+=1

    astrWindowNames.append(strName)

astrClustNames = []
for i in agClustSeqs:
    astrClustNames.append(i.id)

dictCounts = Counter(astrWindowNames)
atupCounts = dictCounts.items()
aCounts = []

for tup in atupCounts:
    aCounts.append(tup[1])

#######################################################
#Tabs on Windows file

sClustNames = set(astrClustNames)
sWindowNames = set(astrWindowNames)

print("VF," + str(len(agWindowSeqs)) + "," + str(iTM) +  "," + str(iQM) + "," + str(len(setHasTM)) + "," + str(len(setHasQM))  +  "," + str(len(sClustNames.difference(sWindowNames))))


print "Total number of windows:"
print len(agWindowSeqs)

print "Breakdown of windows:"
print "(N, Number of genes with N windows)"
for x in Counter(aCounts).items():
    print x

#######################################################
#Check for missing genes




print "Unique genes in clustered file:"
print len(sClustNames)

print "Unique genes in windows file:"
print len(sWindowNames)


print "Number of genes in original clustered file that have 0 windows:"
print len(sClustNames.difference(sWindowNames))

if(args.bList):
    print "List of genes that have 0 windows:"
    for x in sClustNames.difference(sWindowNames):
        print x

