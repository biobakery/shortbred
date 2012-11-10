#!/usr/bin/python

import Bio
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

import re
import sys
import datetime

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
    SeqIO.write( seq, sys.stdout, "fasta")


    
