#!/usr/bin/env python3
#Output phred scores for each index and read file
# Argparse isn't currently on, just put the file name in
#
#
#
#testfile reads.txt.gz


read1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
index1="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
index2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
read2="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

import gzip
#import argparse

#def get_args():
#    parser = argparse.ArgumentParser(description="Get kmer size and file name")
#    parser.add_argument("-c", help="kmer coverage limit", required=True)
#    parser.add_argument("-f", help="input file name", required=True)
#    parser.add_argument("-o", help="output file name", required=True)
#    return parser.parse_args()

#args=get_args()


#INFILE = "{}".format(args.f)
#OUTFILE= "{}".format(args.o)

def convert_phred(letter):
    """Converts a single character into a phred score"""
    ascii = ord('letter')-33
    return(ascii)
file="/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"

def populate_array(file):
    mean_scores = [0.0] * 101

    with gzip.open(file,'rt') as fh:
        i = 0
        LN = 0
        SN= LN/4
        for line in fh:
            LN += 1
            i+=1
            line = line.strip('\n')
            if i%4 == 0:
                #print(line)
                L=0
                #print(LN)
                for letter in line:
                    #print(L)
                    #print(convert_phred(letter))
                    mean_scores[L] += convert_phred(letter)
                    L += 1
    #print(LN/4)
    real_mean_scores=[x/(LN/4) for x in mean_scores]
    #print(real_mean_scores)
    return real_mean_scores
real_mean_scores = populate_array(file)
a = 0
print("# Base Pair    Mean Quality Score")
for a in range(len(real_mean_scores)):
    print(a, real_mean_scores[a], sep="        ")
