#!/usr/bin/env python



import sys
import csv
import argparse
import subprocess
import re

import src.process_blast
pb = src.process_blast

import src.make_windows
mw = src.make_windows

###############################################################################
#ARGUMENTS

parser = argparse.ArgumentParser(description='ShortBRED Identify \n This program produces a set of markers for your genes of interest.')

#INPUT
#goi file
#ref file
#(goi blast results)
#(ref blast results)

parser.add_argument('--goi', type=str, dest='sGOIProts', help='Enter the path and name of the genes of interest file (proteins).')
parser.add_argument('--ref', type=str, dest='sRefProts', help='Enter the path and name of the reference database protein file.')

parser.add_argument('--goiblast', type=str, default = "tmp/goiresults.blast", dest='sGOIBlast', help='Enter the path and name of the blast results from the goi db.')
parser.add_argument('--refblast', type=str, default = "tmp/refresults.blast", dest='sRefBlast', help='Enter the path and name of the blast results from the refrence db.')
parser.add_argument('--markerlength', type=int, default=20, dest='iMLength', help='Enter marker length')

#OUTPUT
#marker file
#map to centroids
parser.add_argument('--markers', type=str, default="markers.faa", dest='sMarkers', help='Enter name and path of marker output file')
parser.add_argument('--cmap', type=str, default="gene-centroid.uc", dest='sMap', help='Enter name and path of centroids-gene file')


#PARAMETERS
parser.add_argument('--threads', type=int, default=1, dest='iThreads', help='Enter the number of threads to use.')
parser.add_argument('--totlength', default = 200, type=int, dest='iTotLength', help='Enter the maximum length for the combined windows for a gene. Default is 200')
parser.add_argument('--id',default = .90, type=float, dest='dID', help='Enter the identity cutoff. Examples: .90, .85, .10,...')
parser.add_argument('--len', default = .10, type=float, dest='dL', help='Enter maximum length for region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
parser.add_argument('--runblast', default = "True", type=str, dest='sRunBlast', help='Set equal to \'False\' to skip Blast, and use existing Blast output')


args = parser.parse_args()

#Clean genes? Check for nuc sequences?

#############################################################################################
#Cluster genes 

log = open("log.txt", "w")
log.write("ShortBRED log \n")




subprocess.check_call(["usearch6", "--cluster_fast", str(args.sGOIProts), "--uc", args.sMap, "--id", ".95","--centroids", "tmp/clust.faa"])


#Remember to output map of genes to their centroids, this is the uc file in usearch, will have to clean it.


#############################################################################################
#BLAST



#print(args.sRunBlast)

if(args.sRunBlast == "True"):

    #Make blastdb's of clustered input genes.
    #MAKE SURE THAT THESE SLASHES DO NOT CAUSE PROBLEMS ON MAC OR WINDOWS.
    subprocess.check_call(["makeblastdb", "-in", "tmp/clust.faa", "-out", "tmp/goidb", "-dbtype", "prot"])
    
    
    #Check the blast DB's
    #subprocess.check_call(["blastdbcheck","-db","tmp/refdb"])
    
    
    
    #Blast input genes against self, and the reference db.
    subprocess.check_call(["blastp", "-query", "tmp/clust.faa", "-db", "tmp/goidb", "-out", "tmp/goiresults.blast", "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])
    
    subprocess.check_call(["makeblastdb", "-in", str(args.sRefProts),"-out", "tmp/refdb/refblastdb", "-dbtype", "prot", "-logfile", "largedb.txt"])    
    
    subprocess.check_call(["blastp", "-query", "tmp/clust.faa", "-db", "tmp/refdb/refblastdb", "-out", "tmp/refresults.blast", "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])

else:
    print "Skipped BLAST."
##################################################################################################
#PROCESS BLAST RESULTS, COUNT OVERLAP BETWEEN GENES (CENTROIDS) AND "HITS"

#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
#dictGOIGenes has form (genename, "AMNLJI....")
#dictRefCounts,dictGOICounts have form (genename,[list of overlap counts for each AA])

dictGOIGenes = pb.getGeneData(open("tmp/clust.faa"))
print "Finding overlap with reference database..."
dictRefCounts = pb.getOverlapCounts(args.sRefBlast, args.dID, 0, args.dL, 0, 0)
print "Finding overlap with goi database..."
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
print "Found True Markers..."
atupQuasiMarkers1 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
print "Found first set of Quasi Markers..."
atupQuasiMarkers2 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
print "Found second set of Quasi Markers..."


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

premarkers = open('tmp/premarkers.txt', 'w')

for key in dictGOIGenes:
    if key in setHasMarkers:
        strGeneName = ">" + key + "_TM"    
        premarkers.write(strGeneName  + '\n')
        premarkers.write(re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)  + '\n')
        iCount = iCount+1
    
premarkers.close()

##################################################################################
#PRINT WINDOWS

dictGeneWindows = mw.getGeneWindows (open('tmp/premarkers.txt'))
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

