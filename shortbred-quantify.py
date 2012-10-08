import sys
import argparse

parser = argparse.ArgumentParser(description='ShortBRED Quantify \n This program takes a set of protein family markers, and produces a relative abundance table.')


parser.add_argument('--markers', type=str, dest='sMarkers', help='Enter the path and name 
of the genes of interest file (proteins).')

parser.add_argument('--wgs', type=str, dest='sGOIProts', help='Enter the path and name of the genes of interest file (proteins).')




#Make a database from the markers

strDBName = sMarkers + ".udb"

subprocess.check_call(["usearch6", "--makeudb_usearch", args.sMarkers, "--output", strName])



#Use usearch to checkj for hits (usearch global)




#Use usearch to checkj for hits (usearch local)
