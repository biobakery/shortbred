#########################################################################
# Jim Kaminski
# Huttenhower Lab
# 2/3/2016
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
    * Convert the first --numspaces (args.iSpace) to "_"\'s \
    * Add an unique ID when two seqs have the exact same name. \
    * Replace characters like "*,[,:" with _ . ',formatter_class=RawTextHelpFormatter) 
parser.add_argument("--numspaces", default = 2, type=int, dest='iSpacesToChange')
args = parser.parse_args()

dictHeaderCounts = {}
reBadChars=re.compile(r'[\\\/\*\=\:\'\[\]\.\;]')

for strLine in sys.stdin:
    strLine = strLine.strip()
    if strLine[0]==">":
        strHeader = strLine[1:].replace(" ","_",args.iSpacesToChange)
        strHeader = re.sub(reBadChars,"_",strHeader)
        dictHeaderCounts[strHeader] = dictHeaderCounts.get(strHeader,0)+1
        if dictHeaderCounts[strHeader] > 1:
            strHeader = strHeader + "___Copy"+str(dictHeaderCounts[strHeader]).zfill(4)
        print ">" + strHeader
    else:
        print strLine
        
        
        
