#!/usr/bin/env python
import argparse
import sys
import os
from Bio import SeqIO as SeqIO

# Given a file with a bunch of IDs, this script extracts a subset of sequences 
# from a FASTA file.

# Ideal usage: $0 -f x.fa -i y.ids -r 0

# Parse Arguments
parser=argparse.ArgumentParser(description="Extract a subset of FASTA \
                                sequences from a FASTA file")
parser.add_argument("fasta", help="The FASTA file from which a subset \
                                is to be extracted")
parser.add_argument("ids", help="IDs of sequences to be extracted")
parser.add_argument("-r","--repeat", help="extraction will now be done \
                    multiple times if an ID is repeated in the IDs list",
                    action="store_true")
args=parser.parse_args()

# Store values from CMD line args to variables
faFile=args.fasta
idFile=args.ids
repeat=args.repeat

# Ensure FASTA and IDs file exist
if(not os.path.exists(faFile)):
    print("FASTA file {0} does not exist, please check input".format(faFile))
    #Graceful exit
    sys.exit(1)

if(not os.path.exists(idFile)):
    print("IDs file {0} does not exist, please check input".format(idFile))
    #Graceful exit
    sys.exit(1)
    
# Read IDs into an array
with open(idFile) as idData:
    idArr=idData.readlines()
# Create a hash with number of occurrences if repeat is enabled
if(repeat):
    idsDict=dict()
    for i in idArr:
        idsDict[i]=idsDict.get(i,0)+1
        
# Bio.SeqIO object reads FASTA file
with open(faFile,"rU") as faHndl:
    # For each FASTA record
    for record in SeqIO.parse(faHndl, "fasta"):
        # If it is part of the required subset
        if(record.id in idArr):
            # And if repeat pickup is enabled
            if(repeat):
                # Output the sequence as many time as the ID
                # is repeated
                for j in range(0,idsDict[record.id]):
                    print(">{0}\n{1}\n\n".format(record.id,record.seq))
            # If repeat is not enabled, just print once
            else:
                print(">{0}\n{1}\n\n".format(record.id,record.seq))
