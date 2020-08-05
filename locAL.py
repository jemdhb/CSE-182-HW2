import argparse
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt 

my_parser = argparse.ArgumentParser(description='local alignment params')
my_parser.add_argument('seq_file', metavar='seq_file',
                       type=str,
                       help='the path to the file w/ the seqs')
my_parser.add_argument('m', metavar='match', type=int,
                       help='match score')

my_parser.add_argument('s', metavar= 'mismatch', type=int,
                       help='mismatch score')

my_parser.add_argument('d', metavar= 'indel', type=int,
                       help='indel score')
my_parser.add_argument('-a', nargs='?',metavar='print',
                       help='print alignment')
args = my_parser.parse_args()
#print(args.a)
###check indel order may have to redo all graphs ahhhhhhhhhhhhhhhh
###fixed on localfast.py
seqFile=open(args.seq_file, 'r')
match=args.m
mismatch=args.s
indel=args.d
print(indel)
localQ1= open("localQ1.txt", "w")
averages= open("average.txt", "w")


#input_path = args.seq_file
#def locAL(str seqFile, int match, int mismatch, int indel):
def matchAtIndex(x,y):
    if x==y:
        return match
    return mismatch
def readFilePairs(file, pairnum):

    seqFile= open(file, "r")
    list=['']*pairnum*2
    with seqFile as infile:
        i=0
        j=0
        for line in infile:
            if line.startswith(">"):
                continue
            else:
                if i%2==1:
                    x=line[0: len(line)-1]
                    list[i]+=x
                    i+=1
                else:
                    y=line[0: len(line)-1]
                    list[i]+= y
                    i+=1
                if i>pairnum*2: break

    return list

def readFile(file):
    x=""
    y=""
    #parse two sequences
    with seqFile as infile:
        i=0
        for line in infile:
            if line.startswith(">"):
                i+=1
                continue
            else:
                if i==1:
                    x+=line
                    x=x[0: len(x)-1]
                else:
                    y+=line
    return (x,y)
    
def make_matrix(sizex, sizey):
    return [[0]*sizey for i in range(sizex)]

def locAL(x,y):
    #print(x)
    #print(y)
    # create a zero-filled matrix
    #print(len(x))
    A= np.zeros((len(x) + 1, len(y) + 1), np.int)
    strings = np.empty((len(x) + 1, len(y) + 1), dtype='string')
    #print(A)
    #print(strings)
    best = 0
    optloc = [0,0]
    # fill in A in the right order
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            score= A[i,j]
        # the local alignment recurrance rule:
            A[i][j], strings[i,j] = max(
            (A[i][j-1] + indel, "i"), 
            (A[i-1][j] + indel, "d"),
            (A[i-1][j-1] + (match if x[i-1] == y[j-1] else mismatch),"m"), (0,"x"))
            #print(A[i,j])
        # track the cell with the largest score
        if A[i][j] >= best:
            best = A[i][j]
            optloc = [i,j]
           # print(optloc)
    # return the opt score and the best location
    print(best)
    best=np.max(A)
    print(A)
    print(best)
    print(optloc)
    i,j = np.unravel_index(A.argmax(), A.shape)
    print(i,j)
    optloc[0]=i
    optloc[1]=j
    seq1=""
    seq2=""
    if args.a!=None:
        while (optloc[0]>0 and optloc[1]>0 ):
            #if match
            xx= optloc[0]
            yy=optloc[1]
            if A[xx,yy]==0:
                break
            if strings[xx][yy]=="m":
                seq1=x[xx-1]+seq1
                seq2=y[yy-1]+seq2
                optloc[0] -= 1
                optloc[1] -= 1
            elif strings[xx][yy]=="i":
                seq1="-"+seq1
                seq2=y[yy-1]+seq2
                optloc[1] -= 1
            elif strings[xx][yy] == "d":
                seq1=x[xx-1]+seq1
                seq2='-'+seq2
                optloc[0] -= 1
        ii=60
        jj=60
        while ii<len(seq1) and jj<len(seq2):
            localQ1.write(seq1[ii-60:ii]+"\n")
            localQ1.write(seq2[jj-60:jj]+"\n\n")
            ii+=60
            jj+=60
        localQ1.write(seq1[ii-60:len(seq1)]+"\n")
        localQ1.write(seq2[jj-60:len(seq2)]+"\n")
        localQ1.write("high score of: "+str(best)+" with length: "+ str(len(seq1))+"\n")
        print(optloc)
    return best, len(seq1)
def plot2():
    l=readFilePairs("randseq.txt",1000)
    i=0
    j=0
    num=[0]*500
    while i+1<len(l)/2 and j < len(l)/2:
        p=locAL(l[i],l[i+1])
        # print(p)
        i+=2
        num[j]=p[1]
        j+=1
    # setting the ranges and no. of intervals 
    range = (min(num), max(num)) 
    bins = 10 
  
    # plotting a histogram 
    plt.hist(num, bins, range, color = 'green', 
        histtype = 'bar', rwidth = 1) 
  
    # x-axis label 
    plt.xlabel('Alignment Length') 
    # frequency label 
    plt.ylabel('Count') 
    # plot title 
    plt.title('Parameter P1 histogram') 
  
    # function to show the plot 
    plt.show() 
#can manually append each thingy and THEN call plot 3
def trackavg():
    l=readFilePairs("randseq.txt",50)
    i=0
    j=0
    num=[0]*25
    avgL=0
    while i+1<len(l)/2 and j < len(l)/2:
        p=locAL(l[i],l[i+1])
        i+=2
        avgL+=p[1]
        j+=1
    # setting the ranges and no. of intervals 
    average= float(avgL)/j
    averages.write(str(average)+"\n")

#for #1 
pair=readFile(seqFile)
x=pair[0]
y=pair[1]
l=locAL(x,y)
