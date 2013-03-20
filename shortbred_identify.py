#!/usr/bin/env python

#Jim Kaminski
#Huttenhower Lab

import sys
import csv
import argparse
import subprocess
import re
import os
import datetime 

import src
import src.process_blast
pb = src.process_blast

import Bio
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

c_strCDHIT	= "cd-hit"
c_strTIME	= "time"

###############################################################################
#COMMAND LINE ARGUMENTS

parser = argparse.ArgumentParser(description='ShortBRED Identify \n This program produces a set of markers for your genes of interest.')


#INPUT Files
#Mode 1: For initial run: goi genes, ref genes, (map_in)
#Mode 2: For initial run with a blastdb: goi genes, ref fb, (map_in)
#Mode 3: For later runs: goi clust, goi blast output, ref blast output, map 

#(map_in) is required for Mode 3, and is an optional input for Modes 1 and 2. If supplied in modes 1 and 2, ShortBRED will map proteins to families according to the pairs in map_in.

#(map_in) is required in Mode 3 because it contains 
 
parser.add_argument('--goi', type=str, dest='sGOIProts',default= "", help='Enter the path and name of the proteins of interest file.')
parser.add_argument('--ref', type=str, dest='sRefProts',default= "", help='Enter the path and name of the file containing reference protein sequences.')

parser.add_argument('--refdb', type=str, dest='dirRefDB', default= "",help='Enter the path and name for  the blastdb of reference proteins.')

parser.add_argument('--goiblast', type=str, default = "", dest='sGOIBlast', help='Enter the path and name of the blast results from the goi-to-goi search.')
parser.add_argument('--refblast', type=str, dest='sRefBlast', default= "", help='Enter the path and name of the blast results from the goi-to-ref search.')
parser.add_argument('--goiclust', type=str, default ="", dest='sClust', help='Enter the path and name of the clustered genes of interest file.')
parser.add_argument('--map_in', type=str, dest='sMapIn',default="", help='Enter the path and name of the two column file connecting proteins to families.')


# DANIELA NOTE: Put these first 7 elements in one list. Easier to check whether the user 
# had supplied incorrect parameters. Determine control flow based on that.

# DB NOTE: Put in a few examples there. Maybe even provide a few sample fasta files.
#	That will help the user out a great deal. They can see how this will work for different
#	Parameters, etc.


#OUTPUT
#Markers, and marker-to-prot_family map
parser.add_argument('--markers', type=str, default="markers.faa", dest='sMarkers', help='Enter name and path for the marker output file')

#DB NOTE - fix the naming system here
parser.add_argument('--map_out', type=str, default="gene-centroid.uc", dest='sMap', help='Enter name and path for the output map file')


#PARAMETERS
#Clustering
parser.add_argument('--clustid',default = .90, type=float, dest='dClustID', help='Enter the identity cutoff for clustering the genes of interest. Examples: .90, .85, .10,...')
parser.add_argument('--qclustid',default = .90, type=float, dest='dQClustID', help='Enter the identity cutoff for clustering the quasi-markers. Examples: .90, .85, .10,...')

#BLAST Search
parser.add_argument('--threads', type=int, default=1, dest='iThreads', help='Enter the number of threads to use.')
parser.add_argument('--id',default = .90, type=float, dest='dID', help='Enter the identity minimum for a short, high-identity region. Examples: .90, .85, .10,...')
parser.add_argument('--len', default = .10, type=float, dest='dL', help='Enter the length maximum for a short, high-identity region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')

#Markers
parser.add_argument('--markerlength', type=int, default=20, dest='iMLength', help='Enter the minimum marker length.')
parser.add_argument('--totlength', default = 200, type=int, dest='iTotLength', help='Enter the maximum length for the combined markers for a gene. Default is 200')
parser.add_argument('--qthresh', type=int, dest='iThresh',default=30, help='Enter a maximum quasi-score.')


#Tmp Directory
parser.add_argument('--tmpdir', default =os.getcwd() +os.sep + "tmp", type=str, dest='sTmp', help='Set directory for temporary output files.')


args = parser.parse_args()

##############################################################################
#Preliminary: Create temporary folder, open log file

dirTmp = args.sTmp

dirTime = src.check_create_dir( dirTmp + os.sep + "time" )
dirQuasi = src.check_create_dir( dirTmp + os.sep + "quasi" )

#DB Note - You could name the log file after the markers. 

log = open(dirTmp + os.sep +"log.txt", "w")
log.write("ShortBRED log \n" + datetime.date.today().ctime() + "\n MARKER PARAMETERS \n")
log.write("ClustID:" + str(args.dClustID) + "\n")
log.write("MarkerLength:" + str(args.iMLength) + "\n")
log.write("QM Thresh:" + str(args.iThresh) + "\n")


iMode = 0

################################################################################
# Step Zero: Choose program mode based on files supplied by the user.

if (args.sGOIProts!="" and args.sRefProts!=""):
	log.write("Mode 1: Building everything..." + "\n")
	iMode = 1
elif (args.sGOIProts!="" and args.dirRefDB!=""):
	log.write("Mode 2: Using user-supplied ref db..." + "\n")
	iMode = 2
elif (args.sClust !="" and args.sGOIBlast!="" and args.sRefBlast!="" and args.sMapIn!=""):
	log.write("Mode 3: Using existing BLAST and Clustering results..." + "\n")
	iMode = 3
else:
	parser.print_help( )
	raise Exception( "Command line arguments incorrect, must provide either:\n" +
		"\t--goi AND --ref, OR\n" +
		"\t--goi AND --refdb, OR\n" +
		"\t--goiclust AND --goiblast AND --refblast AND --map_in" )

#Set default blastresults, clustfile and map from args if not creating the markers
strClustFile = args.sClust
strBlastRef = args.sRefBlast
strBlastSelf = args.sGOIBlast
strMapFile = args.sMapIn

################################################################################
# Step One: Cluster input genes and make into a blast database.
#
# Save centroids to							"tmp/clust/clust.faa"
# Save blastdb of centroids file to			"tmp/clustdb/goidb"
# If prot-to-fam map does not exist, make	"tmp/clust/clust.map"

#Make directories for clustfile and database.
if(iMode==1 or iMode==2):

	dirClust = src.check_create_dir( dirTmp + os.sep + "clust" )
	dirClustDB = src.check_create_dir( dirTmp + os.sep + "clustdb" )
	   
	strClustFile = dirClust + os.sep + "clust.faa"
	strClustDB = dirClustDB + os.sep + "goidb"
	
	#DB Note: clean up strMap
	strMap = os.path.splitext(strClustFile)[0]	
	
	#If user supplied a map, cluster by those families.
	#Else: cluster the entire file and make map file.
	
	if(args.sMapIn!=""):
		"""
		dictFams = {}
		for astrLine in csv.reader( open(args.sMapIn), csv.excel_tab ):
			#dictFams[protein]=[family]		
			dictFams[str(astrLine[1]).strip()]=str(astrLine[0]).strip()
		
		
		#Create subfolders for user-defined families
		dirFams = dirClust + os.sep + "fams"
		strClustPath = dirClust + os.sep + "clust.faa" 
		
		#DB Note - Probably erase this.
		strClutsDB = dirTmp + os.sep + "clustdb" + os.sep + "goi"
		
	
		if not os.path.exists(dirFams):
			os.makedirs(dirFams)
	
		pb.MakeFamilyFastaFiles(dictFams, args.sGOIProts, dirFams)
		pb.ClusterFams(dirClust)
		"""
	else:
		"""
		Old Code to cluster in usearch
		subprocess.check_call(["time", "-o", dirTime + os.sep + "goiclust.time","usearch6", "--cluster_fast", str(args.sGOIProts), "--uc", strMap + ".uc", "--id", str(args.dClustID),"--centroids", strClustFile,"--maxaccepts","0","--maxrejects","0"])
		strMapFile = dirClust + os.sep + "clust.map"
		pb.printMap(strMap+".uc",strMapFile )
		"""
		#Cluster in cdhit
		subprocess.check_call([c_strTIME, "-o",
			  dirTime + os.sep + "goiclust.time",
			  c_strCDHIT, "-i", str(args.sGOIProts),
				"-o", strClustFile, "-d", "0",
				"-c", str(args.dClustID), "-b", "10","-g", "1"])
		strMapFile = dirClust + os.sep + "clust.map"
		pb.GetCDHitMap( strClustFile + ".clstr", strMapFile )
		
		#Create a folder called "clust/fams", will hold a fasta file for each CD-HIT cluster
		dirFams = src.check_create_dir( dirClust + os.sep + "fams" )
		strClutsDB = dirTmp + os.sep + "clustdb" + os.sep + "goi"
	
		#Make a fasta file for each CD-HIT cluster
		pb.MakeFamilyFastaFiles( strMapFile, str(args.sGOIProts), dirFams)
		
		#Call USEARCH to get consensus seq for each cluster,overwrite the CD-HIT cluster file 
		pb.ClusterFams(dirClust, args.dClustID,strClustFile )
	
	#Make database from goi centroids
	subprocess.check_call([c_strTIME, "-o", dirTime + os.sep + "goidb.time",
		"makeblastdb", "-in", strClustFile, "-out", strClustDB,
		"-dbtype", "prot", "-logfile", dirTmp + os.sep + "goidb.log"])

################################################################################
# Step Two: Create reference database, if not supplied by user.
# (refblast and refdb are blank, ref exists)

# Save blastdb of ref file to	 "tmp/refdb/refdb"
if(iMode==1):
	if (args.sRefBlast == "" and args.dirRefDB == "" and args.sRefProts!=""):
		dirRefDB = src.check_create_dir( dirTmp + os.sep + "refdb" )
		strRefDBPath = dirRefDB + os.sep + "refdb"	
		
		#Recent edits to cluster the reference database ahead of time
		"""	  
		subprocess.check_call(["time", "-o",dirTime + os.sep + "refclust.time","usearch6", "--cluster_fast", str(args.sRefProts), "--uc",  "ref.uc", "--id", str(args.dClustID),"--centroids", strClustFile])		
		"""
		
		subprocess.check_call([c_strTIME, "-o", dirTime + os.sep + "refdb.time",
			"makeblastdb", "-in", str(args.sRefProts),"-out", strRefDBPath,
			"-dbtype", "prot", "-logfile", dirTmp + os.sep +  "refdb.log"])

################################################################################
# Step Three: Run Blast Searches
#
#Blast goi centroids against goidb, save results at  "tmp/blastresults/selfblast.txt"
#Blast goi centroids against refdb, save results at  "tmp/blastresults/refblast.txt"

#DB note - rename selfblast.txt goiblast.txt

#If refdb supplied, use that name. 
if(iMode==1 or iMode==2):
	if (args.dirRefDB!=""):
		strRefDBPath = str(args.dirRefDB)
	
	dirBlastResults = src.check_create_dir( dirTmp + os.sep + "blastresults" )
	strBlastRef = dirBlastResults + os.sep + "refblast.txt"
	strBlastSelf = dirBlastResults + os.sep + "selfblast.txt"

	astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30",
		"-ungapped","-comp_based_stats","F","-window_size","0",
		"-xdrop_ungap","1","-evalue","1e-3","-num_alignments","100000",
		"-max_target_seqs", "100000", "-num_descriptions", "100000",
		"-num_threads",str(args.iThreads)]
	
	#Blast clust file against goidb
	subprocess.check_call(["time", "-o", dirTime + os.sep +"goisearch.time",
		"blastp", "-query", strClustFile, "-db", strClustDB,
		"-out", strBlastSelf] + astrBlastParams)
	
	#Blast clust file against refdb
	subprocess.check_call(["time", "-o", dirTime + os.sep +"refsearch.time",
		"blastp", "-query", strClustFile, "-db",strRefDBPath,
		"-out", strBlastRef] + astrBlastParams)

##################################################################################################
#Step Four: Search BLAST output to find regions of overlap.

#make dictGOIgenes:   dict of (genename, seq) for centroids
#make dictRefCounts: dict of (genename, [overlap count for each AA]),

#Get dict of GeneSeqs, then overlap counts from the Ref and GOI blast results
#dictGOIGenes has form (genename, "AMNLJI....")
#dictRefCounts,dictGOICounts have form (genename,[list of overlap counts for each AA])

dictFams = {}

for astrLine in csv.reader(open(strMapFile),delimiter='\t'):
	dictFams[astrLine[1]]=astrLine[0]
		
dictGOIGenes = pb.getGeneData(open(strClustFile))

#Get short, high-identity hits
sys.stderr.write( "Finding overlap with reference database...")
dictRefCounts = pb.getOverlapCounts(strBlastRef, args.dID, 0, args.dL, 0, 0)
sys.stderr.write( "Finding overlap with goi database...")
dictGOICounts = pb.getOverlapCounts(strBlastSelf, args.dID, 0, args.dL, 0, 0)
dictGOICounts = pb.MarkX(dictGOIGenes,dictGOICounts)

#Get medium, high-identity hits - keep edges
dictBigGOICounts = pb.getOverlapCounts(strBlastSelf, args.dID, args.dL +.01, .80, args.iMLength/2, 0)
#Note on 11/1/2012 - for higher precision, could be done for the reference set as well?

#If a gene has 0 valid hits in the ref OR GOI databases, make an array of 0's
#so the program knows that nothing overlapped with the gene.
for dictCounts in (dictRefCounts, dictGOICounts):
	setNotInCounts = set(dictGOIGenes.keys()).difference(set(dictCounts.keys()))
	
	if len(setNotInCounts)>0:
		for sGene in setNotInCounts:
			dictCounts[sGene] = [0]*len(dictGOIGenes[sGene])

#Get dict of counts for (Ref+GOI)
dictAllCounts = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys())) 

for sGene in setRefGOI:
	aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]),dictBigGOICounts.get(sGene,[0]))]
	dictAllCounts[sGene] = aiSum

###########################################################################
#Step Five: Look in AA overlap arrays in dictAllCounts for TM's and QM's

#TM's - look in [AA overlap] for regions of all O's, at least as long as minLength
#QM's - find region if length (minLength) in [AA overlap] with lowest sum([AA overlap])

#DB Note- double-check to make sure this is the same.

setHasMarkers = pb.CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictAllCounts.keys()), dictAllCounts, args.iMLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)
sys.stderr.write( "Found True Markers...")
atupQuasiMarkers1 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength,args.iThresh, args.iTotLength)

sys.stderr.write( "Found first set of Quasi Markers...")
#atupQuasiMarkers2 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iMLength)
#sys.stderr.write( "Found second set of Quasi Markers...")

#Replace AA's with +'s in True Markers
for key in setHasMarkers:
	if dictAllCounts.has_key(key):
		aiWindow = dictAllCounts[key]
		astrGene = list(dictGOIGenes[key])

		for i in range(len(aiWindow)):
			if aiWindow[i] >= 1:
				astrGene[i] = "+"
		dictGOIGenes[key] = "".join(astrGene)

#####################################################################################################
#Step Six: Cluster the Quasi-Markers. Remap the proteins they represent to the centroid marker for each cluster.

atupQM = atupQuasiMarkers1 
atupQM = sorted(atupQM, key=lambda tup: tup[0])

strQuasiFN = dirQuasi+ os.sep + "quasi.faa"
strQuasiClust = dirQuasi+ os.sep + "quasiclust.faa"
strQuasiMap = dirQuasi+ os.sep + "quasi.map"
fQuasi = open(strQuasiFN,'w')
pb.PrintQuasiMarkers(atupQM,fQuasi)
fQuasi.close()

subprocess.check_call(["cd-hit", "-i", strQuasiFN,"-o",strQuasiClust,
	"-d", "0", "-c", str(args.dQClustID), "-b", "10","-g", "1","-aL","1.0"])
		
pb.GetCDHitMap( strQuasiClust+".clstr", strQuasiMap)

dictQuasiClust = {}
for astrLine in csv.reader( open(strQuasiMap), csv.excel_tab ):
			
			mtchMarker = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',astrLine[1])
			strMarker = mtchMarker.group(1)
			mtchFam = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',astrLine[0])
			strFam = mtchFam.group(1)
			
			dictQuasiClust[strMarker] = strFam

for qckey, qcvalue in dictQuasiClust.items():
	for key in dictFams:
		if (dictFams[key] == qckey):
			dictFams[key] = qcvalue

with open(dirTmp + os.sep + "final.map",'w') as fFinalMap:
	for prot, fam in sorted(dictFams.items(), key = lambda(prot, fam): (fam,prot)):
			fFinalMap.write(fam + "\t" + prot + "\n")

#Print AA with overlap area removed to premarkers.txt
strGeneName = ""
iCount = 0
with open(args.sTmp + os.sep + 'premarkers.txt', 'w') as premarkers:
	for key in dictGOIGenes:
		if key in setHasMarkers:
			strGeneName = ">" + key + "_TM"	
			premarkers.write(strGeneName  + '\n')
			premarkers.write(re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)  + '\n')
			iCount = iCount+1

##################################################################################
##Step Seven: Print the TM's, 

#Print the TM's.
#Go through premarkers.txt, find regions satisying user criteria.
"""
Sample gene in premarkers.txt looks like:
>AAA26613_TM
MT+++++LET+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++S++++++
++++++++CLINETEKFLNIWIESNVSF++++++YKSDLLEYKDT+++++++++++++G+++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++QNIVDSVNEWDLNLK
"""	 

fOut = open(args.sMarkers, 'w') 
iMLength = args.iMLength
iTotLength = args.iTotLength

for gene in SeqIO.parse(open(args.sTmp + os.sep + 'premarkers.txt'), "fasta"):	
	iCount = 1
	iRemSeq = iTotLength
	
	mtch = re.search('\+',str(gene.seq))
	if not mtch:
		strMarker = str(gene.seq)
		geneMarker = SeqRecord(Seq(strMarker[0:min(iRemSeq,len(strMarker))]),id = gene.id +"_#" + str(iCount).zfill(2) + '\n', description = "")
		SeqIO.write(geneMarker, fOut,"fasta")
		iRemSeq = iRemSeq - len(geneMarker.seq)

	else:
		for strMarker in (re.split('\++',str(gene.seq))):
			if (iRemSeq>=iMLength and len(strMarker) >= iMLength ):
				geneMarker = SeqRecord(Seq(strMarker[0:min(iRemSeq,len(strMarker))]),id = gene.id +"_#" + str(iCount).zfill(2) + '\n', description = "")
				SeqIO.write(geneMarker, fOut,"fasta")
				iCount+=1
				iRemSeq = iRemSeq - len(geneMarker)

#Print QM's
#Each QM_tuple has the values: (name, window, overlapvalue)

atupQMFinal = []

for tup in atupQM:
	if tup[0] in dictQuasiClust.values():
		atupQMFinal.append(tup)

with open(args.sMarkers, 'a') as fOut:
	pb.PrintQuasiMarkers(atupQMFinal,fOut)
