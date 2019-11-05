#!/usr/bin/env python3

#Demultiplexing
# ```de-multiplex samples to create 48 FASTQ files that contain acceptable
# index pairs (read1 and read2 for 24 different index pairs), two FASTQ files with index-hopped
# reads-pairs, and two FASTQ files undetermined (non-matching or low quality) index-pairs```
#
# #Test files
# read1="TestReads1.gz"
# read2="TestReads2.gz"
# index1="TestIndex1.gz"
# index2="TestIndex2.gz"
# IndexFile="IndexFile"
# CutoffI=25
# CutoffR=30

#argparse input
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('read1', type=str, help='first reads')
parser.add_argument('index1',type=str, help='first indexes')
parser.add_argument('index2',type=str, help='second indexes')
parser.add_argument('read2',type=str, help='second reads')
parser.add_argument('IndexFile',type=str, help='index list file')
parser.add_argument('CutoffI',type=int, help='qscore cutoff for index')
#parser.add_argument('CutoffR', type=int, help='qscore cutoff for reads')
args = parser.parse_args()

def convert_phred(letter):
    """Converts a single character into a phred score"""
    ascii = ord(letter)-33
    return(ascii)

def qfilter(i1, i2):
    """checks if index quality score passes cutoff"""
    #print(i1)
    for letter in i1:
        letter=convert_phred(letter)
        #print(letter)
        if letter < args.CutoffI:
            return(False)
    #print(i2)
    for letter in i2:
        letter=convert_phred(letter)
        #print(letter)
        if letter < args.CutoffI:
            return(False)

def reverse(index):
    """Find reverse compliment of barcode"""
    letters=list(index)
    letters.reverse()
    comp=[]
    for letter in letters:
        if letter == "A":
            comp.append("T")
        elif letter == "T":
            comp.append("A")
        elif letter == "G":
            comp.append("C")
        elif letter == "C":
            comp.append("G")
        else:
            comp.append("N")
    Rcomp=''.join(comp)
    #print(Rcomp)
    return(Rcomp)

##Make dictionary of index pairs
indexpairs={}
with open(args.IndexFile,'rt') as fh:
    LN=0
    for line in fh:
        line = line.strip('\n')
        if LN > 0:
            parts=line.split()
            indexpairs[parts[4]]=reverse(parts[4])
           # print(parts[4])
        LN+=1
        #print("hi")
#print(indexpairs)

import gzip
##open all infiles and most outfiles (matching will open later)
with gzip.open(args.read1,'rt') as R1, \
gzip.open(args.read2, "rt") as R2, \
gzip.open(args.index1, "rt") as I1,\
gzip.open(args.index2, "rt") as I2, \
gzip.open("LowQualityR1.fastq.gz", "wt+") as LR1,\
gzip.open("LowQualityR2.fastq.gz", "wt+") as LR2, \
gzip.open("UndeterminedR1.fastq.gz", "wt+") as UR1,\
gzip.open("UndeterminedR2.fastq.gz", "wt+") as UR2, \
gzip.open("IndexHoppedR1.fastq.gz", "wt+") as IR1, \
gzip.open("IndexHoppedR2.fastq.gz", "wt+") as IR2:
##make counters to track reads in each category
    n=0
    lc=0
    ih=0
    un=0
    ma={}
    ma2=0
##put 4 lines at a time from each file into tuple
    while True:
        header=R1.readline().strip()
        if header=="":
            break
        seq=R1.readline().strip()
        plus=R1.readline().strip()
        qscore=R1.readline().strip()
        R1tup=(header, seq, plus, qscore)
        #print(R1tup)

        header2=R2.readline().strip()
        seq2=R2.readline().strip()
        plus2=R2.readline().strip()
        qscore2=R2.readline().strip()
        R2tup=(header2, seq2, plus2, qscore2)

        header3=I1.readline().strip()
        seq3=I1.readline().strip()
        plus3=I1.readline().strip()
        qscore3=I1.readline().strip()
        I1tup=(header3, seq3, plus3, qscore3)

        header4=I2.readline().strip()
        seq4=I2.readline().strip()
        plus4=I2.readline().strip()
        qscore4=I2.readline().strip()
        I2tup=(header4, seq4, plus4, qscore4)

        n+=1
        if n==1:
            print("1 record read")
##Write lines to appropriate file, add indexes to header
        #print("R1:", R1tup, "R2:", R2tup, "I1:", I1tup, "I2:", I2tup, "\n")
        if qfilter(I1tup[3], I2tup[3]) != False:
            if I1tup[1] in indexpairs:
                if I2tup[1]==indexpairs[I1tup[1]]:
                    #print(I1tup[1], I2tup[1], "match")
                    with gzip.open("Matched_"+I1tup[1]+"-"+I2tup[1]+"_R1.fastq.gz", "at+") as MR1, \
                    gzip.open("Matched_"+I1tup[1]+"-"+I2tup[1]+"_R2.fastq.gz", "at+") as MR2:
                        MR1.writelines([R1tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                        R1tup[1], "\n", R1tup[2], "\n", R1tup[3], "\n"])
                        MR2.writelines([R2tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                        R2tup[1], "\n", R2tup[2], "\n", R2tup[3], "\n"])
                    ma.setdefault(I1tup[1]+"-"+indexpairs[I1tup[1]], 0)
                    ma[I1tup[1]+"-"+indexpairs[I1tup[1]]] +=1
                    ma2+=1
                elif I2tup[1] in indexpairs.values():
                    #print(I1tup[1], I2tup[1], "index hop")
                    IR1.writelines([R1tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                    R1tup[1], "\n", R1tup[2], "\n", R1tup[3], "\n"])
                    IR2.writelines([R2tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                    R2tup[1], "\n", R2tup[2], "\n", R2tup[3], "\n"])
                    ih+=1
                else:
                    #print("unknown")
                    UR1.writelines([R1tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                    R1tup[1], "\n", R1tup[2], "\n", R1tup[3], "\n"])
                    UR2.writelines([R2tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                    R2tup[1], "\n", R2tup[2], "\n", R2tup[3], "\n"])
                    un+=1
            else:
                #print("unknown")
                UR1.writelines([R1tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                R1tup[1], "\n", R1tup[2], "\n", R1tup[3], "\n"])
                UR2.writelines([R2tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
                R2tup[1], "\n", R2tup[2], "\n", R2tup[3], "\n"])
                un+=1
        else:
            #print("low quality")
            LR1.writelines([R1tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
            R1tup[1], "\n", R1tup[2], "\n", R1tup[3], "\n"])
            LR2.writelines([R2tup[0]+(I1tup[1]+"-"+I2tup[1]), "\n", \
            R2tup[1], "\n", R2tup[2], "\n", R2tup[3], "\n"])
            lc+=1

##Output stats
fih=(ih/n)*100
flc=(lc/n)*100
fun=(un/n)*100
fma=(100-fih-flc-fun)
print("Total records:", n)
print("% index hopped:", fih)
print("% low quality:", flc)
print("% undetermined:", fun)
print("% matched:", fma)
print("matched pairs breakdown:")
for x in ma:
    print(x, "|", ((ma[x]/ma2)*100), "% of matched reads", "|", (fma/ma[x]),"% of all records")
