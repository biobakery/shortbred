#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 11:56:09 2012

@author: jim
"""

#!/usr/bin/env python
import csv
import urllib
import xml.dom.minidom
import sys
import re
import random

import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import Entrez

Entrez.email = 'jim.kaminski@gmail.com'
iNum = int(sys.argv[1])
##############################################################
#Read in the data, aSeqs is where all the info will be stored
astrStartSeqs = []
for seq in SeqIO.parse(sys.stdin, "fasta"):
  astrStartSeqs.append(seq.id)

for i in range(len(astrStartSeqs)-1):
    if astrStartSeqs[i]=="1QCA":
        #print "yes"
        del astrStartSeqs[i]

#Pick random sequences
astrSeqs = []
sRandInts = set()

while len(sRandInts) <iNum:
    sRandInts.add(random.randint(0,len(astrStartSeqs)-1))

for x in sRandInts:
    astrSeqs.append([astrStartSeqs[x],"","","",""])







#For each sequence:
#   Get the Protein UID
#   Use the protein ID to get CDS Info
i = 0
for x in astrSeqs:
    #print x
    
    handle = Entrez.esearch(db="protein", term=x[0])
    record = Entrez.read(handle)
    astrIDs = record["IdList"]    
    
    
    handle = Entrez.efetch(db="protein", id=record["IdList"] ,rettype="gb")
    record = SeqIO.read(handle, "genbank")
    
    for y in record.features:
        if y.type == "CDS":
            strRefID = ""
            strSeqStart= ""
            strSeqEnd = ""
            #print record.id
            #print y.qualifiers["coded_by"]
            strCDS = str(y.qualifiers["coded_by"])

            #Pull Info from CDS section
            mtchName = re.search(r'\[\'(.*)\:([0-9]*)\.\.([0-9]*)',strCDS)
            strRefID = mtchName.group(1)
            strSeqStart = mtchName.group(2)
            strSeqEnd = mtchName.group(3)
            
            
            #Knock out complement text
            mtchComp = re.search(r'complement\(',strRefID)
            if mtchComp:
                fComp = True
            else:
                fComp = False
                
            strRefID = strRefID.replace("complement(", "")
            strRefID = strRefID.strip("\"\:")
            
            #print strRefID, strSeqStart, strSeqEnd
            astrSeqs[i][1] = strRefID
            astrSeqs[i][2] = strSeqStart
            astrSeqs[i][3] = strSeqEnd
    i+=1

    """        
    #get nuc id
    handle = Entrez.esearch(db="nuccore", term=x)
    record = Entrez.read(handle)
    """
    #Does it matter if it's the complement strand?
    #print aastrCDSInfo
    
###############################################################
#



i=0
for x in astrSeqs:
    if x[1]!= "":
        handle = Entrez.esearch(db="nuccore", term=x[1])
        record = Entrez.read(handle)
        strID = record["IdList"]
        
        astrSeqs[i][4] = strID
    i +=1
    

aseqNucSeqs = []

#print astrSeqs
i=0
for x in astrSeqs:
    if x[4] != "":
        handle = Entrez.efetch(db="nuccore", id=x[4],rettype="fasta", retmode="text",seq_start=x[2],seq_stop=x[3])
        record = SeqIO.read(handle, "fasta")
        record.id = x[0]
        aseqNucSeqs.append(record)


"""
for x in aastrNucInfo:        
    handle = Entrez.efetch(db="nuccore", id=x[0],rettype="fasta", retmode="text",seq_start=x[1],seq_stop=x[2])
    record = SeqIO.read(handle, "fasta")
    aseqNucSeqs.append(record)
    
print "Number with NUC info"
print len(aseqNucSeqs)

iCount = 0
for i in range(len(aseqNucSeqs)-1):
    aseqNucSeqs[i].id = astrSeqs[i]
"""

SeqIO.write(aseqNucSeqs, sys.stdout, "fasta")     


