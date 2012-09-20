#!/usr/bin/python
import sys

import Bio
from Bio.Seq import Seq
from Bio import SeqIO

aseqGenes = []

len 

for seq in SeqIO.parse(sys.stdin, "fastq"):
    aseqGenes.append(seq)
    
#print len(aseqGenes)
    
#for seq in aseqGenes:
SeqIO.write(aseqGenes, sys.stdout, "fasta")