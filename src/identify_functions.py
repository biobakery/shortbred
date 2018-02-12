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
#        c_strTIME, "-o", dirTime + os.sep + "goiclust.time",
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
#        c_strTIME, "-o", dirTime + os.sep + "goidb.time",
        cmdMakeBlastDB, "-in", fastaInput, "-out", pathDB,
        "-dbtype", "prot", "-logfile", dirTmp + os.sep + strLog])
    

def Run_BLAST_Protein_Search(cmdBLASTP,fastaInput,pathDB,txtBlastOut,iThreads=1):
    #strLog = os.path.basename(pathDB) + ".log"
    astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped",
                       "-comp_based_stats","F","-window_size","0",
                       "-xdrop_ungap","1","-evalue","1e-3",
                       "-max_target_seqs", "1000000",
                       "-num_threads",str(iThreads)]
    subprocess.check_call([
#        "time", "-o", dirTime + os.sep +"goisearch.time",
        cmdBLASTP, "-query", fastaInput, "-db", pathDB,
        "-out", txtBlastOut] + astrBlastParams)



##############################################################################
def Get_Overlap_Counts_From_Search(txtSearchResults, dIDcutoff, iLengthMin, dLengthMax, iOffset, bSaveHitInfo,dictFields):

# Makes a dictionary of form (GeneID, [0,0,0,0,1,1...]), where the number indicates
# the number of times an amino acid overlaps with a region in the blast output.

# PARAMETERS
# txtSearchResults - BLAST-formatted (fmt6) file of hits
# dIDcutoff - Minimum ID for using a hit
# iLengthMin - Lower bound for using a hit in AA's
# dLengthMax - Upper bound for using a hit (AlignmentLength/QueryLength)
# iOffset   - Number of AminoAcids at each end of valid hit which do NOT get +1 to their overlap counts
# bSaveHitInfo - When true, this records the name of the *Target* hit, and the start and end of the hit on the *Query*




#Read in the blast output line by line
#When the program finds a new QueryGene:
    #Add the last gene's aiCounts to dictAAOverlapCounts -- (gene, abWindow)
    #Add the last gene's atupHitInfo to dictOverlapInfo
    #Make a new aiCounts
#When the program finds the same QueryGene:
    #If the region in the blast hit satisfies our length and ID parameters
        #Set the corresponding location in abWindow equal to (count+1)

#Add the last gene when you finish the file

    strCurQuery = ""
    dictAAOverlapCounts = {}
    dictOverlapInfo = {}
    atupHitsInfo = []
    aiCounts =[]
    iGeneCount = 0
    iLine =0



    #sys.stderr.write("The file is " + txtSearchResults + "\n")
    for aLine in csv.reader( open(txtSearchResults, 'rU'), delimiter='\t' ):

        iLine+=1

        strQueryID = aLine[dictFields["QueryID"]]
        strSubId = aLine[dictFields["SubjectID"]]
        dIdentity =float(aLine[dictFields["Identity"]])/100
        iAln= int(aLine[dictFields["Aln"]] )
        iMismatch = int(aLine[dictFields["Mismatches"]])
        iGap = int(aLine[dictFields["Gap"]] )
        iQStart = int(aLine[dictFields["QueryStart"]] )
        iQEnd = int(aLine[dictFields["QueryEnd"]])
        iSubStart = int(aLine[dictFields["SubjectStart"]] )
        iSubEnd= int(aLine[dictFields["SubjectEnd"]] )
        deVal = float(aLine[dictFields["eValue"]] )
        dBit = float(aLine[dictFields["Bit"]])
        iQLength = int(aLine[dictFields["QueryLength"]])

        #When the current queryID changes, add the last one, and switch to the next gene.
        if strQueryID != strCurQuery:
            if iLine>1:
                dictAAOverlapCounts.setdefault(strCurQuery, aiCounts)
                dictOverlapInfo.setdefault(strCurQuery.strip(), atupHitsInfo)
            iGeneCount = iGeneCount+1
            strCurQuery = strQueryID

            atupHitsInfo = []
            aiCounts = []
            for i in range(iQLength):
                aiCounts.append(0)

        dMatchLength = (iAln) / float(iQLength)


        #If it's valid overlap:
        #   and 1 to each appopriate spot in the array
        #   save the name of the overlapping gene, and its start and end points
        #       Note: The latter region will include AA's outside what is marked as overlap, we still want that info.

        if (dIdentity >= dIDcutoff) and (iAln >=iLengthMin) and (dMatchLength <= dLengthMax) and (strQueryID!=strSubId):
            iOLCount = 0
            for i in range(iQStart-1+iOffset, iQEnd-iOffset):
                aiCounts[i]=aiCounts[i]+1
                iOLCount += 1
            if(bSaveHitInfo == True):
                tupInfo = strSubId, iQStart, iQEnd, iSubStart,iSubEnd, iOLCount
                atupHitsInfo.append(tupInfo)



    dictAAOverlapCounts.setdefault(strCurQuery.strip(), aiCounts)
    dictOverlapInfo.setdefault(strCurQuery.strip(), atupHitsInfo)
    iGeneCount = iGeneCount+1

    return dictAAOverlapCounts, dictOverlapInfo
    #tupAAandHits = (dictAAOverlapCounts, dictOverlapInfo)


###########################################################################
def Create_True_Markers(SBFamily,iN):
    
    liValidRegions =[]    
    iCount = 0
    iStart = 0
    bInValidRegion = False

    for i in SBFamily.aiOverlapTotal:
        if i==0 and (bInValidRegion==False):
            iStart = iCount
            bInValidRegion = True
        if i > 0 and (bInValidRegion==True):
            iEnd = iCount
            bInValidRegion = False
            if iEnd-iStart >=iN:
                liValidRegions.append([iStart,iEnd])
        iCount+=1
        
    if i==0:
        iEnd=iCount
        if iEnd-iStart >=iN:
            liValidRegions.append([iStart,iEnd])
            
            
    for aiSlice in liValidRegions:
        marker_slice = slice(aiSlice[0],aiSlice[1])
        #print(SBFamily.strConsensusSeq[marker_slice])
        SBFamily.aSBMarkers = SBFamily.aSBMarkers + [SB.SBMarker(strFamily=SBFamily.strFamily,
                    iCount=len(SBFamily.aSBMarkers)+1,strType="TM",strSeq=SBFamily.strConsensusSeq[marker_slice])]
        SBFamily.iTM = SBFamily.iTM + 1
        
    #print(SBFamily.aSBMarkers)
###############################################################################
def IsInHit( aiSeqQuery, aiSeqTarget):
    if (aiSeqQuery[0] >= aiSeqTarget[0] and aiSeqQuery[1] <= aiSeqTarget[1]):
        bInHit = True
    else:
        bInHit = False
    return bInHit
###############################################################################
def Create_JunctionMarkers_And_QuasiMarkers(SBFamily,strMarkerType,iShortRegion=25,iMarkerLen=25,iXlimit=1):

    # This function takes a set of protein sequences that need JM's, builds
    # up to 1,000 JM's, sorts them by length, and returns the longest 3.

    atupJM = []
    c_iMaxMarkers = 1000

    # Set Window Length
    # Check if X's under limit
    # 
    atupPotentialJM = []
    sys.stderr.write("Processing "+ SBFamily.strFamily +" ...\n")
    iSeqLength = len(SBFamily.strConsensusSeq)
    iStart = 1
    iMarkerCount = 0
    bHitEnd = False
    bFoundRegion = False

    # Check if sequence happens to be shorter than min marker length.
    if (iSeqLength < iMarkerLen):
        bHitEnd = True

        # Load the hits against other consensus sequences.
        #atupHitInfo = dictGOIHits[strGene]

        ########################################################################
        # Alternative: Check Reference Hits as Well.
        # Revisit this in the future. Initial results very slightly weaker, likely
        # due to increased reliance on QM's. Perhaps implement this again,
        # but adjust QM algorithm to penalize consensus overlap more than
        # reference overlap.
        #######################################################################

        atupHitInfo = SBFamily.atupOverlapGOI + SBFamily.atupOverlapRef

        # Make a window as long as iShortRegion. Slide through the sequence.
        # If you find a region with total X's under the limit, move to additional
        # checking, described in the "else" step.

        while (bHitEnd == False and iMarkerCount<c_iMaxMarkers):
            strSeq = SBFamily.strConsensusSeq[iStart-1:(iStart-1)+iShortRegion]

            #Check if hit the end of the sequence.
            if ((iStart+iShortRegion-1) >= iSeqLength):
                bHitEnd = True
            # Check if this seq is over the X limit. If yes, move up one.
            elif(strSeq.count("X")>iXlimit):
                iStart=iStart+1
            else:

                # This region is a candidate for a JM.
                 # Proceed to deep checking, initialize bool vars.
                bOverlapsSomeSeq = False
                bCheckedAllTup = False
                iTupCounter = 0
                #sys.stderr.write("Candidate Seq: "+ strSeq + " " + str(strSeq.count("X")) + str(iXlimit)+"\n" )


                # If there are no overlapping hits at all for the consensus sequence, skip deep checking.
                if(len(atupHitInfo)==0):
                    bCheckedAllTup = True
                    bOverlapsSomeSeq = False

                # Cycle through all of the hits corresponding to this sequence. Check if any completely overlap.
                while (bOverlapsSomeSeq == False and bCheckedAllTup == False and len(atupHitInfo)>0):
                    tupHitInfo = atupHitInfo[iTupCounter]
                    bOverlapsSomeSeq = IsInHit([iStart,iStart+iShortRegion-1],[tupHitInfo[1],tupHitInfo[2]])

                    # We will also check for complete overlap of the region one amino acid ahead (bOverlapAhead),
                    # and one AA behind (bOverlapBehind).
                    if (iStart-1>1):
                        bOverlapBehind = IsInHit([iStart-1,iStart+iShortRegion-2],[tupHitInfo[1],tupHitInfo[2]])
                    else:
                        bOverlapBehind = False
                    if ((iStart+iShortRegion) <= iSeqLength):
                        bOverlapAhead = IsInHit([iStart,iStart+iShortRegion],[tupHitInfo[1],tupHitInfo[2]])
                    else:
                        bOverlapAhead = False

                    # Any overlap results in failure.
                    bOverlapsSomeSeq = (bOverlapsSomeSeq or bOverlapBehind or bOverlapAhead)

                    if(bOverlapsSomeSeq == False):
                        iTupCounter+=1
                    if(iTupCounter>=len(atupHitInfo)):
                        bCheckedAllTup = True

                #####################################################################################
                # Finished checking region, either because we had overlap, or went through all
                #  the hits and had no overlap.
                #####################################################################################

                # Case 1: Region was a success --> try to increase its length.
                if(bOverlapsSomeSeq == False):
                    bFoundRegion = True
                    iEnd = iStart + iShortRegion-1
                    iLength = iEnd - iStart +1
                    #sys.stderr.write("Seq survived tuple check: "+ strSeq +"\n")

                    bHitRightEnd=False
                    bHitLeftEnd=False

                    # Add as many AA's to the right as possible.
                    while (iLength< iMarkerLen and bHitRightEnd==False):
                        iLength = iEnd - iStart +1
                        if(iEnd<iSeqLength and SBFamily.strConsensusSeq[iStart-1:iEnd].count("X")<iXlimit):
                            iEnd+=1
                        else:
                            bHitRightEnd=True

                    # Add as many AA's to the left as possible.
                    while (iLength< iMarkerLen and bHitLeftEnd==False):
                        iLength = iEnd - iStart +1
                        if(iStart>1 and SBFamily.strConsensusSeq[iStart-1:iEnd].count("X")<iXlimit):
                            iStart-=1
                        else:
                            bHitLeftEnd=True

                    #Add the tuple of information for the JM, then move the window up (iMarkerLen) and try for more.


                    tupJM = (SBFamily.strFamily, SBFamily.strConsensusSeq[iStart-1:iEnd],len(SBFamily.strConsensusSeq[iStart-1:iEnd]), iStart,iEnd,"Junction Marker")
                    atupPotentialJM.append(tupJM)
                    iMarkerCount +=1
                    #sys.stderr.write("Original Seq: "+ strSeq +"\n")
                    strSeq = SBFamily.strConsensusSeq[iStart-1:iEnd]
                    #sys.stderr.write("Should be the same seq: "+ strSeq +"\n")
                    iStart = iStart+iMarkerLen
                    bOverlapsSomeSeq = False
                else:
                    # Case 2: Region failed, move onto the next possibility.
                    iStart = iStart+1
                    bOverlapsSomeSeq = False

    #######################################################################
    # Finished checking *all* regions.
    # Sort the JM's by length, return the longest 3. Require these to be
    # at least 10% longer than iShortRegion.
    #######################################################################

    atupSortedJM = sorted(atupPotentialJM, key=lambda tup: tup[2], reverse=True)
    atupSortedJM=atupSortedJM[:3]

    
    # ADD CODE HERE TO CONSTRUCT THE JUNCTION MARKERS AND ADD TO SHORTBRED
    
    for tupJM in atupSortedJM:
        if (tupJM[2]>=(iShortRegion*1.1)):
            sbMarker = SB.SBMarker(strFamily=tupJM[0],iCount=len(SBFamily.aSBMarkers)+1,strType=strMarkerType,strSeq=tupJM[1])
            SBFamily.aSBMarkers = SBFamily.aSBMarkers + [sbMarker]
            
            if(strMarkerType=="JM"):
                SBFamily.iJM += 1
            elif(strMarkerType=="QM"):
                SBFamily.iQM += 1
                
            print("Junction marker program added a marker:")
            
    return

###########################################################################
def MarkX_SBFamilies(dictSBFams, dictOverlap):

    for strName in dictSBFams.keys():
        for i in range(len(dictSBFams[strName].strConsensusSeq)):
            if dictSBFams[strName].strConsensusSeq[i]=="X" or dictSBFams[strName].strConsensusSeq[i]=="x":
                dictOverlap[strName][i] = dictOverlap[strName][i] + 9999999

    return dictOverlap

def Test_IdentifyRevisions():
    fastaInput ="/home/jim/PythonProjects/shortbred/example/input_prots.faa"
    fastaRef = "/home/jim/PythonProjects/shortbred/example/ref_prots.faa"
    dirTmp = "/home/jim/PythonProjects/shortbred/tmp/test"
    cmdMakeBlastDB = "makeblastdb"
    cmdDiamond = ""
    
    strAlignmentProgram = "blast"

    dirClustDB = check_create_dir( dirTmp + os.sep + "clustdb" )
    strClustDB = dirClustDB + os.sep + "goidb"
    strRefDB = dirClustDB + os.sep + "refdb"
    
            
    #Get short, high-identity hits in reference proteins
    c_dictBLASTP_Fields = {"QueryID":0,"SubjectID":1,"Identity":2,"Aln":3,
                    "Mismatches":4,"Gap":5,"QueryStart":6,"QueryEnd":7,
                    "SubjectStart":8,"SubjectEnd":9,"eValue":10,"Bit":11,
                    "QueryLength":12}
    
    
    fastaConsensus, mapClusters = Cluster_FASTA_File_And_Make_Consensus_Sequences(
            cmdCDHIT="/home/jim/shortbred_dependecies/cd-hit/cdhit/cd-hit",
            cmdMUSCLE="/home/jim/shortbred_dependecies/muscle/muscle3.8.31_i86linux64",
            fastaInput=fastaInput,
            dClustID=.85,dConsThresh=.90,dirTmp=dirTmp)
    
    dictSBFamilies = Construct_ShortBRED_Families(fastaConsensus,fastaInput,mapClusters)


    if(strAlignmentProgram=="blast"):
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

        #Update ILength later
        sys.stderr.write( "Finding overlap with reference database...\n")
        dictRefCounts, dictRefHits = Get_Overlap_Counts_From_Search(txtSearchResults=txtBlastRef, dIDcutoff=.90, iLengthMin=8,
                                                         dLengthMax=.15,iOffset=0,bSaveHitInfo=True,
                                                         dictFields=c_dictBLASTP_Fields)
        
        #Get high-identity hits of *all lengths* in GOI database
        sys.stderr.write( "Finding overlap with family consensus database...\n")
        dictGOICounts, dictGOIHits  = Get_Overlap_Counts_From_Search(txtSearchResults=txtBlastSelf, dIDcutoff=.90, iLengthMin=8,
                                                                     dLengthMax=1.0, iOffset=0,bSaveHitInfo=True,dictFields=c_dictBLASTP_Fields)

    #print(dictGOICounts['YP_001068559'])
    dictGOICounts = MarkX_SBFamilies(dictSBFamilies,dictGOICounts)
    #print(dictGOICounts['YP_001068559'])
    
    # Assign dicrionaries to the Families here.
    for strFam in dictSBFamilies.keys():
        dictSBFamilies[strFam].aiOverlapGOI=dictGOICounts.get(strFam,[0]*len(dictSBFamilies[strFam].strConsensusSeq))
        # This is actually a tuple, may cause errors later. Check later.
        dictSBFamilies[strFam].atupOverlapGOI= dictGOIHits.get(strFam,[])
        
        dictSBFamilies[strFam].aiOverlapRef=dictRefCounts.get(strFam,[0]*len(dictSBFamilies[strFam].strConsensusSeq))
        dictSBFamilies[strFam].atupOverlapRef=dictRefHits.get(strFam,[])
        dictSBFamilies[strFam].aiOverlapTotal=[iGOI + iRef for iGOI,iRef in zip(dictSBFamilies[strFam].aiOverlapGOI,dictSBFamilies[strFam].aiOverlapRef)]
        
        #Create True Markers
        Create_True_Markers(dictSBFamilies[strFam],iN=8)
        
        
        # Create a junction marker if no true markers are possible
        if(dictSBFamilies[strFam].iTM==0):
            Create_JunctionMarkers_And_QuasiMarkers(SBFamily=dictSBFamilies[strFam],strMarkerType="JM",iShortRegion=25,iMarkerLen=25,iXlimit=1)
        if(dictSBFamilies[strFam].iTM==0 and dictSBFamilies[strFam].iJM==0):
            Create_JunctionMarkers_And_QuasiMarkers(SBFamily=dictSBFamilies[strFam],strMarkerType="QM",iShortRegion=25,iMarkerLen=25,iXlimit=1)
        
        print(dictSBFamilies[strFam])

    
    # Cluster all of the quasi markers and junction markers, and merge ShortBRED families that get mapped to same cluster.
        
    # Print out final fasta file.




Test_IdentifyRevisions()


"""
Next functions:
    Run_BLAST_Search()
    * Think about how to best organize the overlap process. If the output
    from blast is the same as diamond, can use the same system.
    Run_DIAMOND_Search()
"""