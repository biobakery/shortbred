#########################################################################
# Jim Kaminski
# Huttenhower Lab
# 2/7/2016
#########################################################################



"""
This script makes small changes to an input fasta file to format the sequence
ids for ShortBRED. It will do the following:
    * Convert the first --numspaces (args.iSpacesToChange) to "_"'s
    * Add an unique ID when two seqs have the exact same name.
    * Replace characters like "*,[,:" with _ .

"""
import sys
import argparse
import re
from argparse import RawTextHelpFormatter
    
parser = argparse.ArgumentParser(description='This script makes small changes to an input fasta file to format the sequence \
ids for ShortBRED. It will do the following: \
    * Add a unique ID when two seqs have the exact same name. \
    * Cut fasta files with very long names to 251 characters match (240 chars, plus possible "___Copy000X" suffix for duplicates.) \
    * Replace characters like "*,[,:" with _ . ',formatter_class=RawTextHelpFormatter) 
#parser.add_argument("--numspaces", default = 2, type=int, dest='iSpacesToChange')
args = parser.parse_args()

dictHeaderCounts = {}
reBadChars=re.compile(r'[\\\/\*\=\:\'\[\]\.\;]')
c_iMaxIDLength = 240


for strLine in sys.stdin:
    strLine = strLine.strip()
    if strLine[0]==">":
        #strHeader = strLine[1:].replace(" ","_",args.iSpacesToChange)
        strHeader = strLine[1:].replace(" ","_")
        if len(strHeader) > c_iMaxIDLength:
            sys.stderr.write("Warning: The following header was too long, and was shortened:\n>" + strHeader + "\n" )
            strHeader = strHeader[0:c_iMaxIDLength+1]
        strHeader = re.sub(reBadChars,"_",strHeader)
        dictHeaderCounts[strHeader] = dictHeaderCounts.get(strHeader,0)+1
        if dictHeaderCounts[strHeader] > 1:
            strHeader = strHeader + "___Copy"+str(dictHeaderCounts[strHeader]).zfill(4)
        print ">" + strHeader
    else:
        print strLine
        