#######################################################################################
# This file is provided under the Creative Commons Attribution 3.0 license.
#
# You are free to share, copy, distribute, transmit, or adapt this work
# PROVIDED THAT you attribute the work to the authors listed below.
# For more information, please see the following web page:
# http://creativecommons.org/licenses/by/3.0/
#
# This file is a component of the SflE Scientific workFLow Environment for reproducible 
# research, authored by the Huttenhower lab at the Harvard School of Public Health
# (contact Curtis Huttenhower, chuttenh@hsph.harvard.edu).
#
# If you use this environment, the included scripts, or any related code in your work,
# please let us know, sign up for the SflE user's group (sfle-users@googlegroups.com),
# pass along any issues or feedback, and we'll let you know as soon as a formal citation
# is available.
#######################################################################################



import sfle
import sfleoo
import sys

Import( "*" )

oo = sfleoo.ooSfle(  fileDirOutput = fileDirOutput, fileDirTmp = fileDirTmp)



# Constants
#is usearch working ok with my "0.95", what if it drosp it and does 0.97?
c_iThreads				=  12
c_floatClustID			=  "0.95"

c_floatMaxLength		=  "0.10"
c_floatMinID			=  "0.90"


c_intMinWinLength		=  "20"
c_intTotWinLength		= "200"


c_strRandSeqs			= "50000"


# --Input Files--
c_fileInputGenes		= oo.fin( "VFs.faa" )
c_fileSimData			= oo.fin("s_1_1_export.dat")
c_fileNucSeqs			= oo.fin("VFs.fna")
c_RefGenesPFasta		= oo.fin("all_prots.faa")
#If you have a large reference database, please make a softlink to it in the input directory


#-- Tmp Files --
c_fileCleanGenes		= oo.ftmp("clean.faa" )
c_fileClustGenes        = oo.ftmp("clust.faa" )

c_BlastInputToSelf		= oo.ftmp("blast_clust_to_clust.txt")
c_BlastInputToRef		= oo.ftmp("blast_clust_to_ref.txt")
c_BlastDropToDrop		= oo.ftmp("blast_drop_to_drop.txt")

c_fileMarkerData		= oo.ftmp("markerdata.faa")


c_fileWindows			= oo.ftmp("windows.faa")
c_WinResults			= oo.ftmp("wincheck.txt")

c_fileMaqOutOne			= oo.ftmp("maqsim1.out")
c_fileMaqOutTwo			= oo.ftmp("maqsim2.out")
c_fileSyntheticWGS		= oo.ftmp("wgs.fna")


# Log Files
c_fileCountClean		= oo.ftmp("count_clean.txt" )
c_fileWinCheck			= oo.ftmp("wincheck.txt")

# Programs
c_fileProgCheckSeqs		= oo.fsrc("check_sequences.py")
c_fileProgXRegions 		= oo.fsrc("process_blast_int.py")
c_fileProgMakeWindows 	= oo.fsrc("make_windows.py")
c_fileProgFastq2Fasta	= oo.fsrc("fastq2fasta.py")
c_fileProgCheckWindows	= oo.fsrc("check_windows.py")

#############################################################################################
#Clean input sequences, then cluster. 


oo.pipe(c_fileInputGenes, c_fileCleanGenes, c_fileProgCheckSeqs )
oo.ex( c_fileCleanGenes, c_fileCountClean, "grep", outpipe = True,e=">", c="",verbose = True )
oo.ex( c_fileCleanGenes, c_fileClustGenes, "usearch6", id = c_floatClustID, cluster_fast =  c_fileCleanGenes, 
		centroids = c_fileClustGenes, verbose = True, uc= clustermap.uc)
		
#############################################################################################
#Make blastdb's of clustered input genes.
#Blast input genes against self, and the reference db.

oo.makeblastpdb(c_fileClustGenes)

#The default value for xdrop_ungap is 7. I make it lower to increase the odds of failing.
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



#############################################################################################
# Markers

#X out the conserved small regions.
oo.ex([c_BlastInputToRef,c_BlastInputToSelf,c_fileClustGenes], c_fileMarkerData , c_fileProgXRegions, pipe=True, args =
[("--fasta",c_fileClustGenes),("--goi",c_BlastInputToSelf),("--ref",c_BlastInputToRef),("--winlength",c_intMinWinLength)])

# Make first set of windows. 
oo.ex(c_fileMarkerData, c_fileWindows , c_fileProgMakeWindows,pipe=True, args = [("-WinLength",c_intMinWinLength),("-TotLength",c_intTotWinLength)]) 

# Check windows 
oo.ex([c_fileClustGenes,c_fileWindows], c_WinResults , c_fileProgCheckWindows, pipe=True, args =  [("--cf",c_fileClustGenes)])

Default([c_WinResults])

################################################################################################
#Make Synthetic WGS Data from nuc seqs

oo.ex([c_fileNucSeqs, c_fileSimData],[c_fileMaqOutOne,c_fileMaqOutTwo], "maq", args= ["simulate",("-N",c_strRandSeqs),("",c_fileMaqOutOne),("",c_fileMaqOutTwo),("",c_fileNucSeqs), ("",c_fileSimData)], verbose=True)

oo.ex(c_fileMaqOutOne,c_fileSyntheticWGS, c_fileProgFastq2Fasta, pipe=True)

