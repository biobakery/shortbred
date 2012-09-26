#!/usr/bin/env python


#have the arg parse thing up here


import sys
import csv
import argparse
import subprocess

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


args = parser.parse_args()

#Clean genes? Check for nuc sequences?

#Cluster genes 
subprocess.check_call(["usearch6", "--cluster_fast", str(args.sGOIProts), "--id", ".95","--centroids", "tmp/clust.faa"])

#Remember to output map of genes to their centroids, this is the uc file in usearch, will have to clean it.


#############################################################################################
#Make blastdb's of clustered input genes.
#Blast input genes against self, and the reference db.



#X out the conserved small regions.


# Make first set of windows. 