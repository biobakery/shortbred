#!/usr/bin/env python


#have the arg parse thing up here


import sys
import csv
import argparse
import subprocess

import src.process_blast_int

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
parser.add_argument('--goiblast', type=file, dest='fGOIBlast', help='Enter the path and name of the blast results from the goi db.')
parser.add_argument('--refblast', type=file, dest='fRefBlast', help='Enter the path and name of the blast results from the refrence db.')

#OUTPUT
#marker file
#map to centroids


#PARAMETERS
parser.add_argument('--threads', type=int, default=1, dest='iThreads', help='Enter the number of threads to use.')

args = parser.parse_args()

#Clean genes? Check for nuc sequences?

#Cluster genes 
subprocess.check_call(["usearch6", "--cluster_fast", str(args.sGOIProts), "--id", ".95","--centroids", "tmp/clust.faa"])


#Remember to output map of genes to their centroids, this is the uc file in usearch, will have to clean it.


#############################################################################################
#BLAST

#Make blastdb's of clustered input genes.
#MAKE SURE THAT THESE SLASHES DO NOT CAUSE PROBLEMS ON MAC OR WINDOWS.
subprocess.check_call(["makeblastdb", "-in", "tmp/clust.faa", "-out", "tmp/goidb"])
subprocess.check_call(["makeblastdb", "-in", str(args.sRefProts),"-out", "tmp/refdb"])

#ADD FEATURE FOR MORE PROCESSORS
#Blast input genes against self, and the reference db.
subprocess.check_call(["blastp", "-query", "tmp/clust.faa", "-db", "tmp/goidb", "-out", "tmp/goiresults.blast", "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])

subprocess.check_call(["blastp", "-query", "tmp/clust.faa", "-db", "tmp/refdb", "-out", "tmp/refresults.blast", "-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped","-comp_based_stats","F","-window_size","0", "-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000", "-max_target_seqs", "100000", "-num_descriptions", "100000","-num_threads",str(args.iThreads)])

"""
oo.blastp([c_fileClustGenes,c_fileClustGenes], c_BlastInputToSelf,
			outfmt="6 std qlen",matrix="PAM30",ungapped="",
			comp_based_stats="F",window_size="0", xdrop_ungap=1,
			evalue=1e-3,num_alignments=100000,max_target_seqs=100000,
			num_descriptions=100000,num_threads = 2,makedb=False)


Default([c_BlastInputToSelf])

####################################################################################
#TRY NOT TO EDIT THIS - TAKES SEVERAL HOURS TO RUN
oo.blast([c_fileClustGenes,c_RefGenesPFasta], c_BlastInputToRef,prog = "blastp", 
			outfmt="6 std qlen",matrix="PAM30",ungapped="",
			comp_based_stats="F",window_size="0",
			evalue=1e-3,num_alignments=100000,max_target_seqs=100000,
			num_descriptions=100000,num_threads = c_iThreads,makedb=True)
######################################################################################



#X out the conserved small regions.


# Make first set of windows. 
"""

#######################################################################################################
#PROCESS BLAST RESULTS, COUNT OVERLAP BETWEEN GENES(CENTROIDS) AND "HITS"

#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
#dictGOIGenes has form (genename, "AMNLJI....")
#dictRefCounts,dictGOICounts have form (genename,[list of overlap counts for each AA])
dictGOIGenes = getGeneData(args.fGOIFasta)
dictRefCounts = getOverlapCounts(args.fRefBlast)
dictGOICounts = getOverlapCounts(args.fGOIBlast)

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
dictBoth = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys())) 

for sGene in setRefGOI:
    aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]))]
    dictBoth[sGene] = aiSum

###########################################################################
#CHECK FOR MARKER WINDOWS
#"marker windows": windows of length N that do not overlap with hits in 
#the goi or reference database


setHasMarkers = CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictBoth.keys()), dictBoth, args.iWinLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)

