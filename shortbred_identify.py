#!/usr/bin/env python



import sys
import csv
import argparse
import subprocess
import re
import os

import src.process_blast
pb = src.process_blast

import src.make_windows
mw = src.make_windows

###############################################################################
#ARGUMENTS

parser = argparse.ArgumentParser(description='ShortBRED Identify \n This program produces a set of markers for your genes of interest.')

############################################################################
#INPUT Files
#For initial run: goi genes, ref genes, (famlist)
#For runs after BLAST and Clustering: goiblast, refblast, clust, (famlist)
#(famlist) is an optional specification of protein families
 
parser.add_argument('--goi', type=str, dest='sGOIProts',default= "", help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--refdb', type=str, dest='dirRefDB', default= "",help='Enter the directory where the reference database is stored')
parser.add_argument('--ref', type=str, dest='sRefProts',default= "", help='Enter the path and name of the reference database protein file.')
parser.add_argument('--map_in', type=str, dest='sMapIn',default="", help='Enter the path and name of the two column file connecting proteins to families.')

parser.add_argument('--goiblast', type=str, default = "", dest='sGOIBlast', help='Enter the path and name of the blast results from the goi db.')
parser.add_argument('--refblast', type=str, dest='sRefBlast', default= "", help='Enter the path and name of the blast results from the refrence db.')
parser.add_argument('--goiclust', type=str, default ="", dest='sClust', help='Enter the path and name of clustered file.')

#OUTPUT
#Markers, and marker-to-prot_family map
parser.add_argument('--markers', type=str, default="markers.faa", dest='sMarkers', help='Enter name and path of marker output file')
parser.add_argument('--map_out', type=str, default="gene-centroid.uc", dest='sMap', help='Enter name and path of centroids-gene file')


#PARAMETERS
parser.add_argument('--markerlength', type=int, default=20, dest='iMLength', help='Enter marker length')
parser.add_argument('--threads', type=int, default=1, dest='iThreads', help='Enter the number of threads to use.')
parser.add_argument('--totlength', default = 200, type=int, dest='iTotLength', help='Enter the maximum length for the combined windows for a gene. Default is 200')
parser.add_argument('--id',default = .90, type=float, dest='dID', help='Enter the identity cutoff. Examples: .90, .85, .10,...')
parser.add_argument('--len', default = .10, type=float, dest='dL', help='Enter maximum length for region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
parser.add_argument('--runblast', default = "True", type=str, dest='sRunBlast', help='Set equal to \'False\' to skip Blast, and use existing Blast output')
parser.add_argument('--tmpdir', default =os.getcwd() +os.sep + "tmp", type=str, dest='sTmp', help='Set directory for temporary output files.')


args = parser.parse_args()
##############################################################################

print os.getcwd()

dirTmp = args.sTmp
print dirTmp

#Make out and tmp directories.
if not os.path.exists("out"):
    os.makedirs("out")
if not os.path.exists(dirTmp):
    os.makedirs(dirTmp)


log = open(dirTmp + os.sep + "log.txt", "w")
log.write("ShortBRED log \n")

iMode = 0

################################################################################
# Step Zero: Choose program mode based on args supplied by user
if (args.sRefProts!=""):
    log.write("Mode 1: Building everything..." + "\n")
    iMode = 1
elif (args.sGOIProts!="" and args.dirRefDB!=""):
    log.write("Mode 2: Using user-supplied ref db." + "\n")    
    iMode = 2
elif (args.sClust !="" and args.sGOIBlast!="" and args.sRefBlast!="" and args.sMapIn!=""):
    log.write("Mode 3: Using existing BLAST and Clustering results." + "\n")
    iMode = 3
else:
    print "Command line arguments incorrect."

print "Running in mode",iMode,"..."


################################################################################
# Step One: Cluster input genes and make into a blast database.
#
# Save clustered file to                    "tmp/clust/clust.faa"
# Save blastdb of clustered file to         "tmp/clustdb/goidb"
# If prot-to-fam map does not exist, make   "tmp/clust/clust.map"

#Make directories for clustfile and database.
if(iMode==1 or iMode==2):
    dirClust = dirTmp + os.sep + "clust"
    dirClustDB = dirTmp + os.sep + "clustdb"

    
    if not os.path.exists(dirClust):
        os.makedirs(dirClust)
    if not os.path.exists(dirClustDB):
        os.makedirs(dirClustDB)
     
       
    strClustFile = dirClust + os.sep + "clust.faa"
    strClustDB = dirClustDB + os.sep + "goidb"
    strMap = os.path.splitext(strClustFile)[0]    
    
    
    
    #If user defined protein families, cluster by family.
    #Else: cluster the entire file.
    
    if(args.sMapIn!=""):
        dictFams = {}
        for astrLine in csv.reader( open(args.sMapIn), csv.excel_tab ):
            #dictFams[family]=[protein]        
            dictFams[str(astrLine[1]).strip()]=str(astrLine[0]).strip()
        
        
        #Create subfolders for user-defined families
        dirFams = dirClust + os.sep + "fams"
        strClustPath = dirClust + os.sep + "clust.faa"    
        strClutsDB = dirTmp + os.sep + "clustdb" + os.sep + "goi"
        
    
        if not os.path.exists(dirFams):
            os.makedirs(dirFams)
    
        pb.MakeFamilyFastaFiles(dictFams, args.sGOIProts, dirFams)
        pb.ClusterFams(dirClust)
    else:
        subprocess.check_call(["usearch6", "--cluster_fast", str(args.sGOIProts), "--uc", strMap + ".uc", "--id", ".95","--centroids", strClustFile])
        strMapFile = dirClust + os.sep + "clust.map"
        pb.printMap(strMap+".uc",strMapFile )
    
    #Make goi database
    subprocess.check_call(["makeblastdb", "-in", strClustFile, "-out", strClustDB, "-dbtype", "prot", "-logfile", dirTmp + os.sep + "goidb.log"])


################################################################################
# Step Two: Create reference database, if not supplied by user (refblast and refdb are blank, ref exists)
#
# Save blastdb of ref file to     "tmp/refdb/refdb"
if(iMode==1):
    if (args.sRefBlast == "" and args.dirRefDB == "" and args.sRefProts!=""):
        dirRefDB = dirTmp + os.sep + "refdb"
        if not os.path.exists(dirRefDB):
            os.makedirs(dirRefDB)
        strRefDBPath = dirRefDB + os.sep + "refdb"    
        
        subprocess.check_call(["makeblastdb", "-in", str(args.sRefProts),"-out", strRefDBPath, "-dbtype", "prot", "-logfile", dirTmp + os.sep +  "refdb.log"])

################################################################################
# Step Three: Run Blast Searches
#
#Blast clustered file against goidb, save results at  "tmp/blastresults/selfblast.txt"
#Blast clustered file against refdb, save results at  "tmp/blastresults/refblast.txt"

#If refdb supplied, use that name
if(iMode==1 or iMode==2):
    if (args.dirRefDB!=""):
        strRefDBPath = str(args.dirRefDB)
    
    dirBlastResults = dirTmp + os.sep + "blastresults"
    strBlastRef = dirBlastResults + os.sep + "refblast.txt"
    strBlastSelf = dirBlastResults + os.sep + "selfblast.txt"
    
    
    if not os.path.exists(dirBlastResults):
        os.makedirs(dirBlastResults)
    
    #Blast clust file against goidb
    subprocess.check_call(["blastp", "-query", strClustFile, "-db", strClustDB, "-out", strBlastSelf, "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])
    
    #Blast clust file against refdb
    subprocess.check_call(["blastp", "-query", strClustFile, "-db",strRefDBPath, "-out", strBlastRef, "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])


##################################################################################################
#Part Four: Process BLAST results, find regions of overlap.

#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
#dictGOIGenes has form (genename, "AMNLJI....")
#dictRefCounts,dictGOICounts have form (genename,[list of overlap counts for each AA])

#Get blastresults, clustfile and map from args if only creating the markers
if (args.sClust!=""):
    strClustFile = args.sClust
if (args.sRefBlast!=""):
    strBlastRef = args.sRefBlast
if (args.sGOIBlast!=""):
    strBlastSelf = args.sGOIBlast
if (args.sMapIn!=""):
    strMapFile = args.sMapIn

dictFams = {}
for strLine in csv.reader(open(strMapFile),delimiter='\t'):
    dictFams[strLine[1]]=strLine[0]
        


dictGOIGenes = pb.getGeneData(open(strClustFile))
sys.stderr.write(strBlastRef)
sys.stderr.write(strBlastSelf)
sys.stderr.write( "Finding overlap with reference database...")
dictRefCounts = pb.getOverlapCounts(strBlastRef, args.dID, 0, args.dL, 0, 0)
sys.stderr.write( "Finding overlap with goi database...")
dictGOICounts = pb.getOverlapCounts(strBlastSelf, args.dID, 0, args.dL, 0, 0,dictFams)
dictBigGOICounts = pb.getOverlapCounts(strBlastSelf, args.dID, args.dL +.01, .70, args.iMLength/2, 0)


#If a gene has 0 valid hits in the ref database, make an array of 0's
#so the program knows that nothing overlapped with the gene
setGOINotInRef = set(dictGOIGenes.keys()).difference(set(dictRefCounts.keys()))

if len(setGOINotInRef)>0:
    for sGene in setGOINotInRef:
        dictRefCounts[sGene] = [0]*len(dictGOIGenes[sGene])

#If a gene has 0 valid hits in the GOI database (unlikely), make an 
#array of 0's so the program knows that nothing overlapped with the gene    

setGOINoHits = set(dictGOIGenes.keys()).difference(set(dictGOICounts.keys()))

if len(setGOINoHits)>0:
    for sGene in setGOINoHits:
        dictGOICounts[sGene] = [0]*len(dictGOIGenes[sGene])


#Get dict of counts for (Ref+GOI)
dictAllCounts = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys())) 

for sGene in setRefGOI:
    aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]),dictBigGOICounts.get(sGene,[0]))]
    dictAllCounts[sGene] = aiSum

###########################################################################
#CHECK FOR MARKER WINDOWS, BUILD QUASI-MARKERS FOR LEFTOVER GENES
#"marker windows": windows of length N that do not overlap with hits in 
#the goi or reference database


setHasMarkers = pb.CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictAllCounts.keys()), dictAllCounts, args.iMLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)
sys.stderr.write( "Found True Markers...")
atupQuasiMarkers1 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
sys.stderr.write( "Found first set of Quasi Markers...")
atupQuasiMarkers2 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
sys.stderr.write( "Found second set of Quasi Markers...")


#Replace AA's with X's in True Markers
for key in setHasMarkers:
        if dictAllCounts.has_key(key):
            aiWindow = dictAllCounts[key]
            strGene = list(dictGOIGenes[key])

            for i in range(0,len(aiWindow)):
                if aiWindow[i] >= 1:
                    strGene[i] = "X"
            strGene = "".join(strGene)
            dictGOIGenes[key] =strGene

#Add in the QuasiMarkers
#for key in dictQuasiMarkers:
#    dictGOIGenes[key]=dictQuasiMarkers[key][0]
   
   
#Print out the genes
strGeneName = ""
iCount = 0

premarkers = open(args.sTmp + os.sep + 'premarkers.txt', 'w')

for key in dictGOIGenes:
    if key in setHasMarkers:
        strGeneName = ">" + key + "_TM"    
        premarkers.write(strGeneName  + '\n')
        premarkers.write(re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)  + '\n')
        iCount = iCount+1
    
premarkers.close()

##################################################################################
#PRINT WINDOWS

#Print the True Markers.

#Read in the data from "premarkers.txt",
#It looks like the example below:
"""
>AAA26613_TM
MTXXXXXLETXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXSXXXXXX
XXXXXXXXCLINETEKFLNIWIESNVSFXXXXXXYKSDLLEYKDTXXXXXXXXXXXXXGXXXXXXXXXXXXXXXXXXXXX
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXQNIVDSVNEWDLNLK
"""     
#For each gene,look at the stretches of sequence not X'ed out. 
#If it's long to be a marker, convert it to a seq object and print it out, and
#add the length to iLength. Don't print if iLength would be over the
#total_marker_length limit.

fOut = open(args.sMarkers, 'w') 

for gene in SeqIO.parse(open(args.sTmp + os.sep + 'premarkers.txt'), "fasta"):
    
    iCount = 1
    iLength = 0
    
    for strMarker in (re.split('X*',str(gene.seq))):
        if (len(strMarker)>=args.iMLength and (iLength + len(strMarker)) <= args.iTotLength ):
            geneMarker = SeqRecord(Seq(strMarker),id = ">" + gene.id +"_#" + str(iCount).zfill(2) + '\n', description = "")
            SeqIO.write(geneMarker, fOut,"fasta")
            iCount+=1
            iLength = iLength + len(strMarker)            
            
            




#dictGeneWindows = mw.getGeneWindows (open(args.sTmp + os.sep + 'premarkers.txt'))
#dictSplitWindows = mw.splitGenes(dictGeneWindows, args.iTotLength)
#mw.printWindows(dictSplitWindows, args.sMarkers, args.iMLength, args.iTotLength)
#mw.printQM(dictQuasiMarkers, dictQuasiMarkers2, args.sMarkers)
#print atupQuasiMarkers1
#print atupQuasiMarkers2 

atupQM = atupQuasiMarkers1 + atupQuasiMarkers2
atupQM = sorted(atupQM, key=lambda tup: tup[0])

iCounter = 0
strName = ""


fOut = open(args.sMarkers, 'a')

for tup in atupQM:
    if str(tup[0]) != strName:
        iCounter =1  
    else:
        iCounter+=1
    fOut.write(">" + str(tup[0]) + "_QM" + str(tup[2]) + "_#" +str(iCounter).zfill(2) + '\n')
    fOut.write(str(tup[1]) + '\n')
    strName = str(tup[0])





















"""
if(args.sRunBlast == "True" or args.sRunBlast == "OnlyBlast" ):

   

    
    #Make blastdb's of clustered input genes.
    subprocess.check_call(["makeblastdb", "-in", strClustPath, "-out", strClustDB, "-dbtype", "prot", "-logfile", "smalldb.txt"])
    
    makeblastdb -in /home/jim/tmp/tmp/clust/clust.faa -out /home/jim/tmp/tmp/clustdb/goi

    
    #Blast input genes against self
    subprocess.check_call(["blastp", "-query", strClustName, "-db", args.sTmp + os.sep + "goidb", "-out", args.sGOIBlast, "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])
    
    
    
    subprocess.check_call(["makeblastdb", "-in", str(args.sRefProts),"-out", args.sTmp + os.sep + "refblastdb", "-dbtype", "prot", "-logfile", "largedb.txt"])    
    
    subprocess.check_call(["blastp", "-query", strClustName, "-db", args.sTmp + os.sep + "refblastdb", "-out", args.sRefBlast, "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])
    if(args.sRunBlast == "OnlyBlast"):
        os._exit(1)
else:
    sys.stderr.write( "Skipped BLAST.")
##################################################################################################
#Part Two: PROCESS BLAST RESULTS, COUNT OVERLAP BETWEEN GENES (CENTROIDS) AND "HITS"

#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
#dictGOIGenes has form (genename, "AMNLJI....")
#dictRefCounts,dictGOICounts have form (genename,[list of overlap counts for each AA])
if (args.sClust!=""):
    strClustName = args.sClust



dictGOIGenes = pb.getGeneData(open(strClustName))
sys.stderr.write( "Finding overlap with reference database...")
dictRefCounts = pb.getOverlapCounts(args.sRefBlast, args.dID, 0, args.dL, 0, 0)
sys.stderr.write( "Finding overlap with goi database...")
dictGOICounts = pb.getOverlapCounts(args.sGOIBlast, args.dID, 0, args.dL, 0, 0)
dictBigGOICounts = pb.getOverlapCounts(args.sGOIBlast, args.dID, args.dL +.01, .70, args.iMLength/2, 0)

#If a gene has 0 valid hits in the ref database, make an array of 0's
#so the program knows that nothing overlapped with the gene
setGOINotInRef = set(dictGOIGenes.keys()).difference(set(dictRefCounts.keys()))

if len(setGOINotInRef)>0:
    for sGene in setGOINotInRef:
        dictRefCounts[sGene] = [0]*len(dictGOIGenes[sGene])

#If a gene has 0 valid hits in the GOI database (unlikely), make an 
#array of 0's so the program knows that nothing overlapped with the gene    

setGOINoHits = set(dictGOIGenes.keys()).difference(set(dictGOICounts.keys()))

if len(setGOINoHits)>0:
    for sGene in setGOINoHits:
        dictGOICounts[sGene] = [0]*len(dictGOIGenes[sGene])


#Get dict of counts for (Ref+GOI)
dictAllCounts = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys())) 

for sGene in setRefGOI:
    aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]),dictBigGOICounts.get(sGene,[0]))]
    dictAllCounts[sGene] = aiSum

###########################################################################
#CHECK FOR MARKER WINDOWS, BUILD QUASI-MARKERS FOR LEFTOVER GENES
#"marker windows": windows of length N that do not overlap with hits in 
#the goi or reference database


setHasMarkers = pb.CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictAllCounts.keys()), dictAllCounts, args.iMLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)
sys.stderr.write( "Found True Markers...")
atupQuasiMarkers1 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
sys.stderr.write( "Found first set of Quasi Markers...")
atupQuasiMarkers2 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
sys.stderr.write( "Found second set of Quasi Markers...")


#Replace AA's with X's in True Markers
for key in setHasMarkers:
        if dictAllCounts.has_key(key):
            aiWindow = dictAllCounts[key]
            strGene = list(dictGOIGenes[key])

            for i in range(0,len(aiWindow)):
                if aiWindow[i] >= 1:
                    strGene[i] = "X"
            strGene = "".join(strGene)
            dictGOIGenes[key] =strGene

#Add in the QuasiMarkers
#for key in dictQuasiMarkers:
#    dictGOIGenes[key]=dictQuasiMarkers[key][0]
   
   
#Print out the genes
strGeneName = ""
iCount = 0

premarkers = open(args.sTmp + os.sep + 'premarkers.txt', 'w')

for key in dictGOIGenes:
    if key in setHasMarkers:
        strGeneName = ">" + key + "_TM"    
        premarkers.write(strGeneName  + '\n')
        premarkers.write(re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)  + '\n')
        iCount = iCount+1
    
premarkers.close()

##################################################################################
#PRINT WINDOWS

dictGeneWindows = mw.getGeneWindows (open(args.sTmp + os.sep + 'premarkers.txt'))
dictSplitWindows = mw.splitGenes(dictGeneWindows, args.iTotLength)
mw.printWindows(dictSplitWindows, args.sMarkers, args.iMLength, args.iTotLength)
#mw.printQM(dictQuasiMarkers, dictQuasiMarkers2, args.sMarkers)
#print atupQuasiMarkers1
#print atupQuasiMarkers2 

atupQM = atupQuasiMarkers1 + atupQuasiMarkers2
atupQM = sorted(atupQM, key=lambda tup: tup[0])

iCounter = 0
strName = ""


fOut = open(args.sMarkers, 'a')

for tup in atupQM:
    if str(tup[0]) != strName:
        iCounter =1  
    else:
        iCounter+=1
    fOut.write(">" + str(tup[0]) + "_QM" + str(tup[2]) + "_#" +str(iCounter).zfill(2) + '\n')
    fOut.write(str(tup[1]) + '\n')
    strName = str(tup[0])
    

#concatenate these two lists
#('ABE02101', 'SIFMMGISLTISGLLPQSGF', 20)]
#('ABE02101', 'IGGLLLGLFGNYQKRILLIT', 20)]
"""
