#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 10 19:48:06 2018

@author: jim
"""
import re
import sys
import csv
import argparse
import os
import subprocess
import glob
import math
import shutil

import ShortBRED_Classes as SB

import Bio
from Bio import SeqIO
from Bio import AlignIO
from Bio import Align
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
def check_create_dir( strDir ):

	if not os.path.exists( strDir ):
		os.makedirs( strDir )
	return strDir


def check_file(strPath):
	try:
		with open(strPath):
			pass
	except IOError:
		print( strPath, "does not exist. Please check the path to this file." )
		return

def GetCDHitMap ( fileCluster, strTxtMap):
	setGenes = set()
	strFam = ""
	iLine = 0

	dictFams = {}

	for astrLine in csv.reader(open(fileCluster),delimiter='\t'):
		iLine+=1
		mtch = re.search(r'^>',astrLine[0])
		if (mtch):
			if (len(setGenes)>0 and strFam!=""):
				for strGene in setGenes:
					dictFams[strGene] = strFam
				setGenes=set()
				strFam = ""
		else:
			mtchGene = re.search(r'>(.*)\.\.\.',str(astrLine[1]))
			strGene = mtchGene.groups(1)[0]
			setGenes.add(strGene.strip())
			mtchFam = re.search(r'\.\.\.\s\*',str(astrLine[1]))
			if mtchFam:
				strFam = strGene.strip()

	for strGene in setGenes:
					dictFams[strGene.strip()] = strFam.strip()


	f = open(strTxtMap, 'w')
	for prot, fam in sorted(dictFams.items(), key = lambda a: (a[1], a[0])):
		f.write(fam + "\t" + prot + "\n")
	f.close()

def CheckFastaForBadProtNames(fileFasta):
    reBadChars=re.compile(r'[\\\/\*\=\:\'\[\]\.\;\,]')
    reMarkerTypes = re.compile(r'_([TJQ]M)')
    reWeight = re.compile(r'w=')
    setProtNames = set()
    
    for gene in SeqIO.parse(fileFasta, "fasta"):
        mtchBad = reBadChars.search(gene.id)
        mtchMarkerType = reMarkerTypes.search(gene.id)
        mtchWeight = reWeight.search(gene.id)
        assert (mtchBad == None),("\n\nOne or more of the sequences in your "+
        "input file has an id that ShortBRED cannot use as a valid folder "+
        "name during the clustering step, so ShortBRED has stopped. Please edit  ** "+
        fileFasta + " ** to remove any slashes,asterisks, etc. from the fasta ids. The program utils/AdjustFastaHeadersForShortBRED.py "+
        "in the ShortBRED folder can do this for you.  ShortBRED halted on this gene/protein:" + gene.id)
        assert (gene.id not in setProtNames),("\n\nShortBRED uses the first word of the seq id to identify each "+
        "input sequence, and two or more of your sequences share the same starting word. Please edit ** "+
        fileFasta + " ** to avoid duplicate ids. The program utils/AdjustFastaHeadersForShortBRED.py can add unique identifiers to your data if needed. ShortBRED halted on this gene/protein: " + gene.id)
        assert (mtchMarkerType == None),("\n\nShortBRED uses the strings '_TM','_QM',and '_JM' in the headers for markers, "+
        "and their presence in the original protein ID may cause issues downstream. Please edit ** "+
        fileFasta + " ** . ShortBRED halted on this gene/protein: " + gene.id)
        assert (mtchWeight == None),("\n\nShortBRED uses the string '_w=' in the headers for markers, "+
        "and their presence in the original protein ID may cause issues downstream. Please edit ** "+
        fileFasta + " ** . ShortBRED halted on this gene/protein: " + gene.id)
        
        setProtNames.add(gene.id)


def MakeFamilyFastaFiles( strMapFile, fileFasta, dirOut):
    #Makes a fasta file containing genes for each family in dictFams.
    dictFams = {}
    for strLine in csv.reader(open(strMapFile),delimiter='\t'):
        dictFams[strLine[1]]=strLine[0]
    dictgeneFamilies = {}

    for gene in Bio.SeqIO.parse(fileFasta, "fasta"):
        strFamily = dictFams.get(gene.id,"unassigned")
  

        if strFamily == "unassigned":
            
            # REPLACE WITH LOGGER LATER
            sys.stderr.write("Warning: The sequence " + gene.id + " was not assigned\
 to any cluster by cd-hit, and will not have a marker. This often happens when a\
 sequence is very short (say <8 AA's).\n")
        else:
            aGene = dictgeneFamilies.get(strFamily,list())
            aGene.append(gene)
            dictgeneFamilies[strFamily] = aGene

    for key in dictgeneFamilies.keys():
        strFamFile = dirOut + os.sep + key + ".faa"
        f = open(strFamFile, 'w')
        Bio.SeqIO.write(dictgeneFamilies[key], strFamFile, "fasta")
        f.close()


def ClusterFams(dirClust, dCLustID, strOutputFile, dThresh,strMUSCLE):
#Clusters all of the family files made by MakeFamilyFastaFiles.

    dirFams = dirClust + os.sep + "fams"
    dirCentroids = dirFams+os.sep+"centroids"
    dirUC = dirFams+ os.sep + "uc"

    if not os.path.exists(dirClust):
        os.makedirs(dirClust)
    if not os.path.exists(dirFams):
        os.makedirs(dirFams)
    if not os.path.exists(dirCentroids):
        os.makedirs(dirCentroids)

    if not os.path.exists(dirUC):
        os.makedirs(dirUC)


	#sys.stderr.write( dirCentroids + "\n")
	#sys.stderr.write( str(glob.glob(dirFams+os.sep+'*.faa')) + "\n")
    for fileFasta in glob.glob(dirFams+os.sep+'*.faa'):
        #sys.stderr.write("The file is " + fileFasta + " \n")
        fileClust = dirCentroids + os.sep + os.path.basename(fileFasta)
        fileAlign = dirFams + os.sep + os.path.basename(fileFasta)+".aln"
        strSeqID = os.path.basename(fileFasta)
        strSeqID = strSeqID.replace(".faa", "")

        iSeqCount = 0
        #Count seqs, if more than one, then align them
        for seq in SeqIO.parse(fileFasta, "fasta"):
            iSeqCount+=1
            
        if iSeqCount>1:
            #Call muscle to produce an alignment
            subprocess.check_call([strMUSCLE, "-in", str(fileFasta), "-out", str(fileAlign)])

            # Use BioPython's "dumb consensus" feature to get consensus sequence
            algnFasta = AlignIO.read(str(fileAlign), "fasta")

            seqConsensus =str(AlignInfo.SummaryInfo(algnFasta).dumb_consensus(threshold=dThresh, ambiguous='X'))
            seqConsensus = SeqRecord(Seq(seqConsensus),id=strSeqID)

            SeqIO.write(seqConsensus, str(fileClust), "fasta")


            """
			# We previously used EMBOSS-CONS to produce consensus sequences
			# Call cons or em_cons from the EMBOSS package to produce a consensus sequence
			subprocess.check_call(["cons", "-seq", str(fileAlign), "-outseq", str(fileClust)])
            """
        else:
            shutil.copyfile(fileFasta,fileClust)



    ageneAllGenes = []

    for fileFasta in glob.glob(dirCentroids+os.sep+'*.faa'):
        for gene in SeqIO.parse(fileFasta, "fasta"):
            gene.id = os.path.basename(os.path.splitext(fileFasta)[0])
            ageneAllGenes.append(gene)

    """
	for gene in ageneAllGenes:
		mtch = re.search(r'centroid=(.*)',gene.id)
		if mtch:
			gene.id = mtch.group(1)
		else:
			gene.id = os.path.splitext()
    """
    SeqIO.write(ageneAllGenes, strOutputFile, "fasta")
    return strOutputFile


def Cluster_FASTA_File_And_Make_Consensus_Sequences(cmdCDHIT,cmdMUSCLE,fastaInput,dClustID,dConsThresh,
                                              dirTmp):
    """ Calls CD-Hit to cluster the sequences. Returns 'map' file which
    has two columns [FamilyID,SequenceID]. """
    
    
    dirClust = check_create_dir( dirTmp + os.sep + "clust" )
    dirClustDB = check_create_dir( dirTmp + os.sep + "clustdb" )
    strClustFile = dirClust + os.sep + "clust.faa"
    strClustDB = dirClustDB + os.sep + "goidb"

    sys.stderr.write( "Clustering proteins of interest...\n")
    check_file(fastaInput)

    subprocess.check_call([
#		c_strTIME, "-o", dirTime + os.sep + "goiclust.time",
		  cmdCDHIT, "-i", str(fastaInput),
			"-o", strClustFile, "-d", "0",
			"-c", str(dClustID), "-b", "10","-g", "1"])
    strMapFile = dirClust + os.sep + "clust.map"
    GetCDHitMap( strClustFile + ".clstr", strMapFile )
    sys.stderr.write( "Protein sequences clustered.")

    sys.stderr.write( "Creating folders for each protein family...\n")
    #Create a folder called "clust/fams", will hold a fasta file for each CD-HIT cluster
    dirFams = check_create_dir( dirClust + os.sep + "fams" )
    strClustDB = dirTmp + os.sep + "clustdb" + os.sep + "goi"

    sys.stderr.write( "Making a fasta file for each protein family...\n")
    #Make a fasta file for each CD-HIT cluster
    MakeFamilyFastaFiles( strMapFile, fastaInput, dirFams)
    #Cluster these and make consensus sequences
    fastaConsensus = ClusterFams(dirClust,dClustID,strClustFile,dConsThresh,cmdMUSCLE )
    return([fastaConsensus,strMapFile])


def Construct_ShortBRED_Families(fastaConsensus,fastaInput,mapClust):
    dictSBFamilies = {}
    dictInputSeqs = {}
    dictClust = {}
    
    with open(mapClust,'r') as fileClust:
        for strLine in fileClust:
            strFam, strSeq = strLine.split("\t")[0], strLine.split("\t")[1]
            dictClust[strFam] = list(dictClust.get(strFam,[])) + [strSeq.strip()]
    print(dictClust)
    
    for seq in SeqIO.parse(fastaInput, "fasta"):
        dictInputSeqs[seq.id] = seq.seq
    
    for consensus_seq in SeqIO.parse(fastaConsensus, "fasta"):
        dictSBFamilies[consensus_seq.id] = SB.SBFamily(strFamily=consensus_seq.id,strConsensusSeq=consensus_seq.seq)
        dictMemberSeqs = {}
        for strMemberSeq in dictClust[consensus_seq.id]:
            dictMemberSeqs[consensus_seq.id] = dictInputSeqs[strMemberSeq]
        dictSBFamilies[consensus_seq.id].dictMemberSeqs = dictMemberSeqs
    return dictSBFamilies

def Create_BLAST_Database(cmdMakeBlastDB,fastaInput,pathDB,dirTmp):

    strLog = os.path.basename(pathDB) + ".log"
    sys.stderr.write("Making BLAST database for the family consensus sequences...\n")
    #Make database from goi centroids
    subprocess.check_call([
#		c_strTIME, "-o", dirTime + os.sep + "goidb.time",
        cmdMakeBlastDB, "-in", fastaInput, "-out", pathDB,
        "-dbtype", "prot", "-logfile", dirTmp + os.sep + strLog])
    

def Run_BLAST_Protein_Search(cmdBLASTP,fastaInput,pathDB,txtBlastOut,iThreads=1):
    strLog = os.path.basename(pathDB) + ".log"
    astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped",
                       "-comp_based_stats","F","-window_size","0",
                       "-xdrop_ungap","1","-evalue","1e-3",
                       "-max_target_seqs", "1000000",
                       "-num_threads",str(iThreads)]
    subprocess.check_call([
#		"time", "-o", dirTime + os.sep +"goisearch.time",
        cmdBLASTP, "-query", fastaInput, "-db", pathDB,
        "-out", txtBlastOut] + astrBlastParams)



def Test_ClusterFASTAFileAndMakeFamilies():
    fastaInput ="/home/jim/PythonProjects/shortbred/example/input_prots.faa"
    fastaRef = "/home/jim/PythonProjects/shortbred/example/ref_prots.faa"
    dirTmp = "/home/jim/PythonProjects/shortbred/tmp/test"
    cmdMakeBlastDB = "makeblastdb"

    dirClustDB = check_create_dir( dirTmp + os.sep + "clustdb" )
    strClustDB = dirClustDB + os.sep + "goidb"
    strRefDB = dirClustDB + os.sep + "refdb"
    
    
    fastaConsensus, mapClusters = Cluster_FASTA_File_And_Make_Consensus_Sequences(
            cmdCDHIT="/home/jim/shortbred_dependecies/cd-hit/cdhit/cd-hit",
            cmdMUSCLE="/home/jim/shortbred_dependecies/muscle/muscle3.8.31_i86linux64",
            fastaInput=fastaInput,
            dClustID=.85,dConsThresh=.90,dirTmp=dirTmp)
    
    dictSBFamilies = Construct_ShortBRED_Families(fastaConsensus,fastaInput,mapClusters)
    #print(dictSBFamilies)
    #print(dictSBFamilies["P13246"].dictMemberSeqs)
    Create_BLAST_Database(cmdMakeBlastDB,fastaConsensus, strClustDB,dirTmp)
    Create_BLAST_Database(cmdMakeBlastDB,fastaRef, strRefDB,dirTmp)
    
    
    # Blast searches
    dirBlastResults = check_create_dir( dirTmp + os.sep + "blastresults" )
    txtBlastRef = dirBlastResults + os.sep + "refblast.txt"
    txtBlastSelf = dirBlastResults + os.sep + "selfblast.txt"
    
    Run_BLAST_Protein_Search(cmdBLASTP="blastp",fastaInput=fastaConsensus,
                             pathDB=strClustDB,txtBlastOut=txtBlastSelf,iThreads=1)

    Run_BLAST_Protein_Search(cmdBLASTP="blastp",fastaInput=fastaConsensus,
                             pathDB=strRefDB,txtBlastOut=txtBlastRef,iThreads=1)

Test_ClusterFASTAFileAndMakeFamilies()


"""
Next functions:
    Run_BLAST_Search()
    * Think about how to best organize the overlap process. If the output
    from blast is the same as diamond, can use the same system.
    Run_DIAMOND_Search()
"""