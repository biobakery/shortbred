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
#from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

import re
import sys
import datetime
import argparse

parser = argparse.ArgumentParser(description='This program checks for short seqs, and nucleotide seqs.')
parser.add_argument('--minlength', type=int, dest='iN', help='Enter the the minimum length for a sequence.',default =1)
args = parser.parse_args()

def getSeqs(fileFasta):

    aSeqs = []

    for seq in SeqIO.parse(fileFasta, "fasta"):
        aSeqs.append(seq)

    return aSeqs


aSeqs = getSeqs(sys.stdin)


#fix to see if there are consensus sequences
#delete nucleotide sequences
for x in aSeqs:
    mtch = re.search('[^ACGTN]',str(x.seq))
    if mtch == None:
        aSeqs.remove(x)

#print new list to fasta file
for seq in aSeqs:
	if len(seq) > args.iN:
		SeqIO.write( seq, sys.stdout, "fasta")



