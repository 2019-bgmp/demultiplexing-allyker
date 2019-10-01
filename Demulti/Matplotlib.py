#!/usr/bin/env python3
#NEED CONDA ON
#conda install matplotlib
#change file for each plot
#Index1
#Index2
#Reads1
#Reads2

file="Reads1"
Xarray=[]
Yarray=[]


#print("h")
with open(file, "r") as fh:
    #LN=0
    for line in fh:
        line = line.strip('\n')
        #if LN > 0:
        parts=line.split()
        Xarray.append(float(parts[0]))
        Yarray.append(float(parts[1]))
        #print("hi")
        #LN+=1
#print(Yarray)


import matplotlib.pyplot as mp
mp.scatter(Xarray,Yarray,label='Average Q score, Read2',color='r')
mp.xlabel('Base Pair #')
mp.ylabel('Average Q score')
mp.title('Average Q score at Each Base')
mp.legend()
mp.show()
#print("Hi")
