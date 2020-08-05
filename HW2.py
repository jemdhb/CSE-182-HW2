import argparse
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt 
#argument parser
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
seqFile=open(args.seq_file, 'r')
match=args.m
mismatch=args.s
indel=args.d
#opening files ahead of time so I can write data onto the same file 
#with different functions
localQ1= open("localQ1.txt", "w")
averages= open("average.txt", "w")
#uses > as a marker between two sequences
def readFilePairs(file, pairnum):

    seqFile= open(file, "r")
    #list should be twice as long as the number of pairs
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
                #here for files that have stats at the bottom
                #after all of the sequences
                if i>pairnum*2: break

    return list
#used to read a file with only one pair
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
    
#aligns  sequences locally 
def locAL(x,y):
    # create a zero-filled matrix
    A= np.zeros((len(x) + 1, len(y) + 1), np.int)
    strings = np.empty((len(x) + 1, len(y) + 1), dtype='string')
    #tracks best score
    best = 0
    #tracks location of the best score
    optloc = [0,0]
    # fill in A with scores and strings with directionality 
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
        # should mark the max score in the areas adjacent to 
        # the cell as the next place to go:
            A[i][j], strings[i,j] = max(
            (A[i][j-1] + indel, "i"), 
            (A[i-1][j] + indel, "d"),
            (A[i-1][j-1] + (match if x[i-1] == y[j-1] else mismatch),"m"), (0,"x"))
        # track the cell with the largest score
        if A[i][j] >= best:
            best = A[i][j]
            optloc = [i,j]
    # return the opt score and the best location
    #was having issues with one off indexes with the above tracking
    #method
    best=np.max(A)
    i,j = np.unravel_index(A.argmax(), A.shape)
    optloc[0]=i
    optloc[1]=j
    #used for final backtracking alignments 
    seq1=""
    seq2=""
    #backtracking 
    if args.a!=None:
        while (optloc[0]>0 and optloc[1]>0 ):
            xx= optloc[0]
            yy=optloc[1]
            #exit as soon as a zero score is reached
            if A[xx,yy]==0:
                break
            #if a match or mismatch move diagonally 
            if strings[xx][yy]=="m":
                seq1=x[xx-1]+seq1
                seq2=y[yy-1]+seq2
                optloc[0] -= 1
                optloc[1] -= 1
            #if an indel move vertically 
            elif strings[xx][yy]=="i":
                seq1="-"+seq1
                seq2=y[yy-1]+seq2
                optloc[1] -= 1
            #if a deletion move horozontally
            elif strings[xx][yy] == "d":
                seq1=x[xx-1]+seq1
                seq2='-'+seq2
                optloc[0] -= 1
        #done to make the sequence readable in the file
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
    #return best score and length of best alignment 
    return best, len(seq1)
#used to plot the average length of alignments 
def plot2():
    #read in randomly generated data 
    l=readFilePairs("randseq.txt",1000)
    i=0
    j=0
    num=[0]*500
    #align every randomly generated pair 
    while i+1<len(l)/2 and j < len(l)/2:
        p=locAL(l[i],l[i+1])
        i+=2
        num[j]=p[1]
        j+=1
    # setting the ranges and no. of intervals 
    range = (min(num), max(num)) 
    bins = 10 
  
    # plotting the histogram 
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
#writes averages to a file to be read by another function
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

#for #2 have to call random dnapy first for it to work
plot2()
