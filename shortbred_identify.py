#!/usr/bin/env python
#####################################################################################
#Copyright (C) <2013> Jim Kaminski and the Huttenhower Lab
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of
#this software and associated documentation files (the "Software"), to deal in the
#Software without restriction, including without limitation the rights to use, copy,
#modify, merge, publish, distribute, sublicense, and/or sell copies of the Software,
#and to permit persons to whom the Software is furnished to do so, subject to
#the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies
#or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
#INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A
#PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
#HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
#OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
#SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#
# This file is a component of ShortBRED (Short, Better REad Database)
# authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Jim Kaminski, jjk451@mail.harvard.edu).
#####################################################################################

#Jim Kaminski
#Huttenhower Lab


import sys
import csv
import argparse
from argparse import RawTextHelpFormatter
import subprocess
import re
import os
import time
import datetime
import math



try:
    import Bio
except ImportError:
    print "\nShortBRED was unable to load Biopython. Please check to make sure you have Biopython installed (http://biopython.org/wiki/Main_Page), and that its directory is in your PYTHONPATH. \n"
    sys.exit(1)

import src
import src.process_blast
pb = src.process_blast

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable
from Bio import SeqIO

VERSION="0.9.5"

###############################################################################
#COMMAND LINE ARGUMENTS

parser = argparse.ArgumentParser(description='ShortBRED Identify: \n This program produces a set of markers for your proteins of interest. \n \
The minimum input files required to run the program are: \n \
\t [--goi] 1) A fasta file of proteins, for which you want to build markers. \n \
\t [--ref] 2) A fasta file of reference proteins  \n \
The program will output a file fasta file of markers [--markers]. \n \n\
Example: \n \
\t$ ./ python shortbred_identify.py --goi example/input_prots.faa --ref example/ref_prots.faa --tmp tmp_example \n \
Please note that a large number of comparisons can take several hours to run. [Example: 2,000 prots of interest against a reference set of 5,000,000 prots]  \n \n \
It is also possible to use an existing BLAST protein database as your reference set, or modify parameters from an earlier ShortBRED run. \
Please see the documentation for more details.',formatter_class=RawTextHelpFormatter)
parser.add_argument("--version", action="version", version="%(prog)s v"+VERSION)

#INPUT Files
#Mode 1: For initial run: goi genes, ref genes, (map_in)
#Mode 2: For initial run with a blastdb: goi genes, ref fb, (map_in)
#Mode 3: For later runs: goi clust, goi blast output, ref blast output, map

#(map_in) is required for Mode 3, and is an optional input for Modes 1 and 2. If supplied in modes 1 and 2, ShortBRED will map proteins to families according to the pairs in map_in.
#(map_in) is required in Mode 3 because it contains

grpInput = parser.add_argument_group('Input')
grpInput.add_argument('--goi', type=str, dest='sGOIProts',default= "", help='Enter the path and name of the proteins of interest file.')
grpInput.add_argument('--ref', type=str, dest='sRefProts',default= "", help='Enter the path and name of the file containing reference protein sequences.')
grpInput.add_argument('--refdb', type=str, dest='dirRefDB', default= "",help='Can be specified in place of reference proteins [--ref]. Enter the path and name for a blastdb of reference proteins.')
grpInput.add_argument('--goiblast', type=str, default = "", dest='sGOIBlast', help='Used when modifying existing ShortBRED-Identify results. Enter the path and name of the blast results from the goi-to-goi search.')
grpInput.add_argument('--refblast', type=str, dest='sRefBlast', default= "", help='Used when modifying existing ShortBRED-Identify results. Enter the path and name of the blast results from the goi-to-ref search.')
grpInput.add_argument('--goiclust', type=str, default ="", dest='sClust', help='Used when modifying existing ShortBRED-Identify results. Enter the path and name of the clustered genes of interest file.')
grpInput.add_argument('--map_in', type=str, dest='sMapIn',default="", help='Used when modifying existing ShortBRED-Identify results. Enter the path and name of the two column file connecting proteins to families.')



#OUTPUT
#Markers, and marker-to-prot_family map
grpOutput = parser.add_argument_group('Output')
grpOutput.add_argument('--markers', type=str, default="markers.faa", dest='sMarkers', help='Enter name and path for the marker output file')



grpOutput.add_argument('--map_out', type=str, default="gene-centroid.uc", dest='sMap', help='Enter name and path for the output map file')


#PARAMETERS
grpParam = parser.add_argument_group('Parameters')
# Clustering
grpParam.add_argument('--clustid',default = .85, type=float, dest='dClustID', help='Enter the identity cutoff for clustering the genes of interest. Examples: .90, .85, .10,...')
grpParam.add_argument('--qclustid',default = .90, type=float, dest='dQClustID', help='Enter the identity cutoff for clustering the quasi-markers. Examples: .90, .85, .10,...')
grpParam.add_argument('--consthresh',default = .95, type=float, dest='dConsThresh', help='Enter the consensus threshold for assigning AA\'s in the family alignments to the consensus sequences. The default is .70. Examples: .60, .70, .80,...')

# BLAST Search
grpParam.add_argument('--threads', type=int, default=1, dest='iThreads', help='Enter the number of threads to use.')
grpParam.add_argument('--id',default = .90, type=float, dest='dID', help='Enter the identity minimum for a short, high-identity region. Examples: .90, .85, .10,...')
grpParam.add_argument('--len', default = .15, type=float, dest='dL', help='Enter the length maximum for a short, high-identity region. l=(length hit region)/(length query gene) Examples: .30, .20, .10,... ')
grpParam.add_argument('--minAln', default = 0, type=int, dest='iLenMin', help='Enter the minimum for a short, high-identity region. Examples: 10, 20, 30,... ')

# Markers
grpParam.add_argument('--markerlength', type=int, default=8, dest='iMLength', help='Enter the minimum marker length.')
grpParam.add_argument('--totlength', default = 300, type=int, dest='iTotLength', help='Enter the maximum length for the combined markers for a gene. Default is 200')
grpParam.add_argument('--qthresh', type=int, dest='iThresh',default=1, help='Enter a maximum quasi-score.')
grpParam.add_argument('--qmlength', type=int, dest='iQMlength',default=33, help='Enter a minimum length for QM\'s.')
grpParam.add_argument('--xlimit', dest='iXlimit', default=1,help='Enter the number of Xs to allow in JMs')

# Tmp Directory
grpParam.add_argument('--tmpdir', default ="", type=str, dest='sTmp', help='Set directory for temporary output files.')

# Programs
grpPrograms = parser.add_argument_group('Programs')
grpPrograms.add_argument('--usearch', default ="usearch", type=str, dest='strUSEARCH', help='Provide the path to usearch. Default call will be \"usearch\".')
grpPrograms.add_argument('--muscle', default ="muscle", type=str, dest='strMUSCLE', help='Provide the path to muscle. Default call will be \"muscle\".')
grpPrograms.add_argument('--cdhit', default ="cd-hit", type=str, dest='strCDHIT', help='Provide the path to usearch. Default call will be \"cd-hit\".')
grpPrograms.add_argument('--blastp', default ="blastp", type=str, dest='strBLASTP', help='Provide the path to blastp. Default call will be \"blastp\".')
grpPrograms.add_argument('--makeblastdb', default ="makeblastdb", type=str, dest='strMAKEBLASTDB', help='Provide the path to  makeblastdb. Default call will be to \"blastp\".')

args = parser.parse_args()


##############################################################################
#Preliminary: Check for args, appropriate dependencies 
# Create temporary folder, and open log file.

# Check for args.
if len(sys.argv)==1:
    parser.print_help()
    sys.stderr.write("\nNo arguments were supplied to ShortBRED. Please see the usage information above to determine what to pass to the program.\n")
    sys.exit(1)

# Check dependencies

    #iReturnCode = oCmd.returncode
    #return iReturnCode

print "Checking dependencies..."
src.CheckDependency(args.strUSEARCH,"","usearch")    
src.CheckDependency(args.strBLASTP,"-h","blastp")
src.CheckDependency(args.strMUSCLE,"-h","muscle")
src.CheckDependency(args.strCDHIT,"-h","cdhit")
src.CheckDependency(args.strMAKEBLASTDB,"-h","makeblastdb")

print "Checking to make sure that installed version of usearch can make databases..."
pCmd = subprocess.Popen([args.strUSEARCH,"-help","makeudb_usearch"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
pCmd.communicate()[0]
if (pCmd.returncode==0):
    sys.stderr.write("Usearch appears to be working.\n\n")
else:
    sys.stderr.write("Usearch cannot make a required database. Please check to make sure you have v6.0.307 of usearch or later installed.\n\n")

c_strTIME	= "time"
strCDHIT = args.strCDHIT
strBLASTP = args.strBLASTP
strMUSCLE = args.strMUSCLE




if args.iLenMin==0:
	#If no minimum alignment length is given, use 80% of minimum marker length and round up.
	iLenMin = math.ceil(args.iMLength*(.80))
else:
    iLenMin = args.iMLength




dirTmp = args.sTmp
if(dirTmp==""):
	# dirTmp gets a pid and timestamp. (This is to avoid overwriting files if
	# someone launches multiple instances of the program.)
    dirTmp = ("tmp" + str(os.getpid()) + '%.0f' % round((time.time()*1000), 1))

dirTime = src.check_create_dir( dirTmp + os.sep + "time" )
dirQuasi = src.check_create_dir( dirTmp + os.sep + "quasi" )

#In the future, will use the logging module in Python.
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
	pb.CheckFastaForBadProtNames(args.sGOIProts)    
	dirClust = src.check_create_dir( dirTmp + os.sep + "clust" )
	dirClustDB = src.check_create_dir( dirTmp + os.sep + "clustdb" )

	strClustFile = dirClust + os.sep + "clust.faa"
	strClustDB = dirClustDB + os.sep + "goidb"

	#DB Note: clean up strMap
	strMap = os.path.splitext(strClustFile)[0]

	sys.stderr.write( "Clustering proteins of interest...\n")
	src.check_file(str(args.sGOIProts))
	#Cluster in cdhit
	subprocess.check_call([
#		c_strTIME, "-o", dirTime + os.sep + "goiclust.time",
		  strCDHIT, "-i", str(args.sGOIProts),
			"-o", strClustFile, "-d", "0",
			"-c", str(args.dClustID), "-b", "10","-g", "1"])
	strMapFile = dirClust + os.sep + "clust.map"
	pb.GetCDHitMap( strClustFile + ".clstr", strMapFile )
	sys.stderr.write( "Protein sequences clustered.")

	sys.stderr.write( "Creating folders for each protein family...\n")
	#Create a folder called "clust/fams", will hold a fasta file for each CD-HIT cluster
	dirFams = src.check_create_dir( dirClust + os.sep + "fams" )
	strClutsDB = dirTmp + os.sep + "clustdb" + os.sep + "goi"

	sys.stderr.write( "Making a fasta file for each protein family...\n")
	#Make a fasta file for each CD-HIT cluster
	pb.MakeFamilyFastaFiles( strMapFile, str(args.sGOIProts), dirFams, log)

	sys.stderr.write( "Aligning sequences in each family, producing consensus sequences...\n")
	#Call MUSCLE + EMBOSS_CONS OR DUMB CONSENSUS to get consensus seq for each cluster,overwrite the CD-HIT cluster file
	pb.ClusterFams(dirClust, args.dClustID,strClustFile,args.dConsThresh,args.strMUSCLE )

	sys.stderr.write("Making BLAST database for the family consensus sequences...\n")
	#Make database from goi centroids
	subprocess.check_call([
#		c_strTIME, "-o", dirTime + os.sep + "goidb.time",
		args.strMAKEBLASTDB, "-in", strClustFile, "-out", strClustDB,
		"-dbtype", "prot", "-logfile", dirTmp + os.sep + "goidb.log"])

################################################################################
# Step Two: Create reference database, if not supplied by user.
# (refblast and refdb are blank, ref exists)

# Save blastdb of ref file to	 "tmp/refdb/refdb"
if(iMode==1):
	if (args.sRefBlast == "" and args.dirRefDB == "" and args.sRefProts!=""):
		dirRefDB = src.check_create_dir( dirTmp + os.sep + "refdb" )
		strRefDBPath = dirRefDB + os.sep + "refdb"

		sys.stderr.write("Making BLAST database for the reference protein sequences...\n")
		src.check_file(str(args.sRefProts))
		subprocess.check_call([
#			c_strTIME, "-o", dirTime + os.sep + "refdb.time",
			args.strMAKEBLASTDB, "-in", str(args.sRefProts),"-out", strRefDBPath,
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

	astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30", "-ungapped",
		"-comp_based_stats","F","-window_size","0",
		"-xdrop_ungap","1","-evalue","1e-3",
		"-max_target_seqs", "1000000",
		"-num_threads",str(args.iThreads)]


	#Blast clust file against goidb
	sys.stderr.write( "BLASTing the consensus family sequences against themselves...\n")
	subprocess.check_call([
#		"time", "-o", dirTime + os.sep +"goisearch.time",
		strBLASTP, "-query", strClustFile, "-db", strClustDB,
		"-out", strBlastSelf] + astrBlastParams)

	#Blast clust file against refdb
	sys.stderr.write("BLASTing the consensus family sequences against the reference protein sequences...\n")
	subprocess.check_call([
#		"time", "-o", dirTime + os.sep +"refsearch.time",
		strBLASTP, "-query", strClustFile, "-db",strRefDBPath,
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
log.write("Initial Families:\t" + str(len(dictGOIGenes)) + "\n")

#Get short, high-identity hits in reference proteins
sys.stderr.write( "Finding overlap with reference database...\n")
dictRefCounts, dictRefHits = pb.getOverlapCounts(strBlastRef, args.dID, iLenMin, args.dL, 0,bSaveHitInfo=True)

#Get high-identity hits of *all lengths* in GOI database
sys.stderr.write( "Finding overlap with family consensus database...\n")
dictGOICounts, dictGOIHits  = pb.getOverlapCounts(strBlastSelf, args.dID, iLenMin, 1.0, 0,True)
dictGOICounts = pb.MarkX(dictGOIGenes,dictGOICounts)


"""
#Get medium, high-identity hits - keep edges
dictBigGOICounts, dictBigGOIHits = pb.getOverlapCounts(strBlastSelf, args.dID, args.dL +.01, .80, args.iMLength/2, 0, bSaveHitInfo=True)
#Note on 11/1/2012 - for higher precision, could be done for the reference set as well?
"""

#If a gene has 0 valid hits in the ref OR GOI databases, make an array of 0's
#(length of the gene) so the program knows that nothing overlapped with the gene.

for dictCounts in (dictRefCounts, dictGOICounts):
	setNotInCounts = set(dictGOIGenes.keys()).difference(set(dictCounts.keys()))

	if len(setNotInCounts)>0:
		for sGene in setNotInCounts:
			dictCounts[sGene] = [0]*len(dictGOIGenes[sGene])

#Get dict of counts for (Ref+GOI)
dictAllCounts = {}
setRefGOI = set(dictGOICounts.keys()).union(set(dictRefCounts.keys()))

for sGene in setRefGOI:
	#Old, removed dictforGOICounts
	#aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]),dictBigGOICounts.get(sGene,[0]))]
	#New
	aiSum =[sum(aiCounts) for aiCounts in zip(dictGOICounts.get(sGene,[0]),dictRefCounts.get(sGene,[0]))]
	dictAllCounts[sGene] = aiSum



###########################################################################
#Step Five: Look in AA overlap arrays in dictAllCounts for TM's and QM's

#TM's - look in [AA overlap] for regions of all O's, at least as long as minLength
#QM's - find region if length (minLength) in [AA overlap] with lowest sum([AA overlap])

#DB Note- double-check to make sure this is the same.

setHasMarkers = pb.CheckForMarkers(set(dictGOIGenes.keys()).intersection(dictAllCounts.keys()), dictAllCounts, args.iMLength)
setLeftover = set(dictGOIGenes.keys()).difference(setHasMarkers)
setAll = setHasMarkers.union(setLeftover)

log.write("Families with True Markers:\t" + str(len(setHasMarkers)) + "\n")

sys.stderr.write( "Found True Markers...\n")

#iShort = math.floor(args.iMLength*(.95))
#sys.stderr.write("The Short region is " + str(int(iShort)) )

# Get Junction Markers.
# This line of code should perhaps be moved under "if len(setLeftover)>0:"
atupQuasiMarkers1 = pb.FindJMMarker(setLeftover, dictGOIGenes, dictGOIHits,dictRefHits,iShortRegion = int(math.floor(args.iQMlength*.40)),iXlimit=int(args.iXlimit),iMarkerLen=args.iQMlength)

# This checks to see if any families need a JM or QM.
if len(setLeftover)>0:
	bHasQuasi = True

else:
	bHasQuasi = False
	sys.stderr.write( "No Quasi Markers needed...\n")

if bHasQuasi:
	sys.stderr.write( "Found "+str(len(atupQuasiMarkers1)) +" JM Markers...\n")
	if len(atupQuasiMarkers1) >0:
		for a in zip(*atupQuasiMarkers1):
			setGotQM = a
			break
		setLeftover = setLeftover.difference(setGotQM)

	# Change these lines to determine how QM's are made
	atupQuasiMarkers2 = pb.CheckForQuasiMarkers(setLeftover, dictAllCounts, dictGOIGenes,args.iQMlength,args.iThresh, args.iTotLength)
	atupQuasiMarkers = atupQuasiMarkers1 + atupQuasiMarkers2
	sys.stderr.write("Found " +str(len(atupQuasiMarkers2)) + " QM-Minimals.\n")



#Replace AA's with +'s in True Markers
for key in setHasMarkers:
	if key in dictAllCounts:
		aiWindow = dictAllCounts[key]
		astrGene = list(dictGOIGenes[key])

		for i in range(len(aiWindow)):
			if aiWindow[i] >= 1:
				astrGene[i] = "+"
		dictGOIGenes[key] = "".join(astrGene)

#####################################################################################################
#Step Six: Cluster the Quasi-Markers and Junction Markers. Remap the proteins they represent to the centroid marker for each cluster.

if(bHasQuasi):
	atupQM = atupQuasiMarkers
	atupQM = sorted(atupQM, key=lambda tup: tup[0])

	strQuasiFN = dirQuasi+ os.sep + "quasi.faa"
	strQuasiClust = dirQuasi+ os.sep + "quasiclust.faa"
	strQuasiMap = dirQuasi+ os.sep + "quasi.map"
	fQuasi = open(strQuasiFN,'w')
	pb.PrintQuasiMarkers(atupQM,fQuasi,bDetailed=False,bInitial=True)
	fQuasi.close()

	subprocess.check_call([args.strCDHIT, "-i", strQuasiFN,"-o",strQuasiClust,
		"-d", "0", "-c", str(args.dQClustID), "-b", "8","-g", "1","-aL","1.0"])

	pb.GetCDHitMap( strQuasiClust+".clstr", strQuasiMap)

	dictQuasiClust = {}
	for astrLine in csv.reader( open(strQuasiMap), csv.excel_tab ):

				mtchMarker = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',astrLine[1])
				strMarker = mtchMarker.group(1)
				mtchFam = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',astrLine[0])
				strFam = mtchFam.group(1)

				dictQuasiClust[strMarker] = strFam
	log.write("QM Families, before clustering:\t" +str(len(set(dictQuasiClust.keys()))) + "\n")
	for qckey, qcvalue in dictQuasiClust.items():
		for key in dictFams:
			if (dictFams[key] == qckey):
				dictFams[key] = qcvalue

	with open(dirTmp + os.sep + "final.map",'w') as fFinalMap:
		for prot, fam in sorted(dictFams.items(), key = lambda a: (a[1],a[0])):
				fFinalMap.write(fam + "\t" + prot + "\n")

	log.write("QM Families, after clustering:\t" +str(len(set(dictQuasiClust.values()))) + "\n")

#Print AA with overlap area removed to premarkers.txt
strGeneName = ""
iCount = 0
with open(dirTmp + os.sep + 'premarkers.txt', 'w') as premarkers:
	for key in dictGOIGenes:
		if key in setHasMarkers:
			strGeneName = ">" + key + "_TM"
			premarkers.write(strGeneName  + '\n')
			premarkers.write(re.sub("(.{80})","\\1\n",dictGOIGenes[key],re.DOTALL)  + '\n')
			iCount = iCount+1

##################################################################################
#Step Seven: Print the TM's to a temp file, print information on JM's and QM's.

#Print the TM's.
#Go through premarkers.txt, find regions satisying user criteria.
"""
Sample gene in premarkers.txt looks like:
>AAA26613_TM
MT+++++LET+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++S++++++
++++++++CLINETEKFLNIWIESNVSF++++++YKSDLLEYKDT+++++++++++++G+++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++QNIVDSVNEWDLNLK
"""


dirFrameCheck = src.check_create_dir( dirTmp + os.sep + "framecheck" )
strTmpMarkers = dirFrameCheck + os.sep + "FirstMarkers.faa"
fOut = open(strTmpMarkers, 'w')
iMLength = args.iMLength
iTotLength = args.iTotLength

for gene in SeqIO.parse(open(dirTmp + os.sep + 'premarkers.txt'), "fasta"):
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

################################################################################
# Debugging - Print out information on QM's
if(bHasQuasi):
	# Reload Gene Sequences into dictionary
	dictGOIGenes = pb.getGeneData(open(strClustFile))

	# Output info to txt file
	strQMOut = dirTmp+os.sep+"QMInformation.txt"

	atupQM = pb.UpdateQMHeader(atupQM,dictGOIHits,dictRefHits, strQMOut,dictGOIGenes)



	###############################################################################
	atupQMFinal = []
	# Gather the QM's that survived clustering.
	for tup in atupQM:
		if tup[0] in dictQuasiClust.values():
			atupQMFinal.append(tup)


	# Print the final set of Quasi Markers
	with open(strTmpMarkers, 'a') as fOut:
		pb.PrintQuasiMarkers(atupQMFinal,fOut,True,False)

sys.stderr.write( "\nTmp markers saved to " + strTmpMarkers + "\n")


################################################################################
# Step Eight: Blastx nuc versions of markers against the GOI db (prots)
# Process the results, flag any marker that hits something not in its family.


"""
# Reverse translate each candidate marker (back_translate)

strFrameNucs = dirFrameCheck + os.sep + "MarkersAsNucs.fna"
subprocess.check_call(["backtranseq", "-sequence", strTmpMarkers, "-out", strFrameNucs])


I will change the above function to be a Python function so that we do not need
EMBOSS. "back_translate" was deprecated from BioPython, so it was not an option.
Perhaps I can find an old version online, and copy the function in here and cite
the authors. If not, I can make something simple using the BioPython functions.


#*** BLASTX the whole shebang against family consensus sequences


if (iMode ==3):
    # Make a db for the clustered goi file if needed.
	strClustFile = args.sClust
	strClustDB = dirTmp + os.sep + "clustdb" + os.sep + "goidb"

	subprocess.check_call(["makeblastdb", "-in", strClustFile, "-out", strClustDB,
			"-dbtype", "prot", "-logfile", dirTmp + os.sep + "goidb.log"])

# Blastx crashes when attempting ungapped alignments on short reads
if args.iMLength <30:
	astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30",
	 #"-ungapped", "-xdrop_ungap","1",
	 "-evalue","1e-3",
		"-max_target_seqs", "100000",
		"-num_threads",str(args.iThreads)]
else:
    astrBlastParams = ["-outfmt", "6 std qlen", "-matrix", "PAM30",
	"-ungapped",
		"-xdrop_ungap","1","-evalue","1e-3",
		"-max_target_seqs", "100000",
		"-num_threads",str(args.iThreads)]

#New code for trying out new version of blastx. Throw this out later.
astrBlastParams = ["-outfmt", "6 std qlen"]


#Should I include something with "frame_shift_penalty"?


strBlastNucsToGOI = dirFrameCheck + os.sep + "NucsToGOI.blast"
#Blast clust file against goidb
sys.stderr.write( "BLASTing the backtranslated marker nucleotides against family consensus sequences (blastx)...\n")
subprocess.check_call(["blastx", "-query", strFrameNucs, "-db", strClustDB,
	"-out", strBlastNucsToGOI] + astrBlastParams)

strOffTargetHits = dirFrameCheck + os.sep + "OffTargetHits.txt"
setProblemMarkers = pb.CheckOutOfFrame (strBlastNucsToGOI,.95, 32,dictFams, strOffTargetHits)

log.write("Number flagged as potential off-target markers: " + str(len(setProblemMarkers)) + "\n")

strOffTargetMarkers = dirFrameCheck + os.sep + "ProblemMarkerList.txt"

with open(strOffTargetMarkers, 'w') as fOut:
	for strMarkerName in setProblemMarkers:
		fOut.write(strMarkerName + "\n")

# Remove any marker for which any rev-trans has a significant off-target hit

# Track marker statistics (#TM's,etc.) as processing, and print to log.
"""
iTM = 0
iQM = 0
iQMJunction = 0

setMarkerFamilies = set()
setAllProtFamilies = set(dictFams.values())

#TEMP CHANGE
#Adding this line to keep off-target markers
setProblemMarkers = set()

with open(args.sMarkers,'w') as fOut:
	for gene in SeqIO.parse(open(strTmpMarkers), "fasta"):
		if gene.id not in setProblemMarkers:
			mtchProtStub = re.search(r'(.*)_(.M)[0-9]*_\#([0-9]*)',gene.id)
			strMarkerProt = mtchProtStub.group(1)
			if (gene.id.find("_TM")>0):
				iTM+=1
			if (gene.id.find("_QM")>0):
				iQM+=1
			if (gene.id.find("_JM")>0):
				iQMJunction+=1
			setMarkerFamilies.add(strMarkerProt)
			SeqIO.write(gene, fOut,"fasta")

sys.stderr.write( "\nProcessing complete! Final markers saved to " + args.sMarkers + "\n")
sys.stderr.write( "\nNOTE: Please open and read the log file before using the markers\
. Warnings about individual sequences will appear in the log.\n\n")

iQMMinimal = iQM
iMarkers = iTM + iQMJunction + iQMMinimal

log.write("Total Markers:\t" + str(iMarkers) + "\n")
log.write("True Markers:\t " + str(iTM) + "\n")
log.write("QM-Junctions:\t " + str(iQMJunction) + "\n")
log.write("QM-Minimals:\t " + str(iQMMinimal) + "\n\n")

log.write("Families with Markers:\t " + str(len(setMarkerFamilies)) + "\n")
log.write("Families without Markers:\t " + str(len(setAllProtFamilies.difference(setMarkerFamilies))) + "\n")


