import argparse
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt 
#arg parser
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
localQ1= open("localFast.txt", "w")
#reads files with 2 sequences
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
#space efficent version of local that only 
# finds the maxScore and index
def findMax(x,y):
    #two lines that we will track
    currLine= np.zeros(len(x)+1)
    prevLine=np.zeros(len(x)+1)
    best = 0
    optloc = [0,0]
    #move line by line, without the matrix
    for i in range(1, len(x)+1):
        #swap line values 
        prevLine=currLine
        currLine= np.zeros(len(x)+1)
        for j in range(1, len(y)+1):
            #grab max value from nearby
            currLine[j] = max(
            (currLine[j-1] + indel), 
            (prevLine[j] + indel),
            (prevLine[j-1] + (match if x[i-1] == y[j-1] else mismatch)), (0))
            #when a better value is found, update optloc and best
            if currLine[j] >= best:
                best = currLine[j]
                optloc = [i,j]
    return [optloc, best]
#run on size reduced strings in order to be space efficent  
def locAL(x,y):
    #tracks scores
    A= np.zeros((len(x) + 1, len(y) + 1), np.int)
    #tracks movement
    strings = np.empty((len(x) + 1, len(y) + 1), dtype='string')
    best = 0
    optloc = [0,0]
    # fill in A with scores
    for i in range(1, len(x)+1):
        for j in range(1, len(y)+1):
            score= A[i,j]
        # put max score and direction to match its location:
            A[i][j], strings[i,j] = max(
            (A[i][j-1] + indel, "i"), 
            (A[i-1][j] + indel, "d"),
            (A[i-1][j-1] + (match if x[i-1] == y[j-1] else mismatch),"m"), (0,"x"))
        # track the cell with the largest score
        if A[i][j] >= best:
            best = A[i][j]
            optloc = [i,j]
           # print(optloc)
    # return the opt score and the best location
    best=np.max(A)
    #optimal start should be the bottom right corner of the matrix
    optloc[0]=len(x)
    optloc[1]=len(y)
    seq1=""
    seq2=""
    #backtrack
    if args.a==None:
        #align as normal, since this is a linear subspace
        while (optloc[0]>0 and optloc[1]>0 ):
            #if match
            xx= optloc[0]
            yy=optloc[1]
            if A[xx,yy]==0:
                break
            #match/mismatch= diagonal movement 
            if strings[xx][yy]=="m":
                seq1=x[xx-1]+seq1
                seq2=y[yy-1]+seq2
                optloc[0] -= 1
                optloc[1] -= 1
            #indel = vertical movement 
            elif strings[xx][yy]=="i":
                seq1="-"+seq1
                seq2=y[yy-1]+seq2
                optloc[1] -= 1
            #deletion = horizontal movement
            elif strings[xx][yy] == "d":
                seq1=x[xx-1]+seq1
                seq2='-'+seq2
                optloc[0] -= 1
        #outputs code to file in a readable manner
        ii=60
        jj=60
        while ii<len(seq1) and jj<len(seq2):
            localQ1.write(seq1[ii-60:ii]+"\n")
            localQ1.write(seq2[jj-60:jj]+"\n\n")
            ii+=60
            jj+=60
        localQ1.write(seq1[ii:len(seq1)]+"\n")
        localQ1.write(seq2[jj:len(seq2)])
        print(optloc)
    return best, len(seq1)
#reads in sequences
p= readFile(seqFile)
#inputs sequences into findmax
endIndices=findMax(p[0],p[1])
#reduce sequences to end at findmax result 
# and then reverse the sequence
rev1= p[0]
rev1= rev1[0:endIndices[0][0]+1]
rev2= p[1]
rev2= rev2[0:endIndices[0][1]+1]
#call find max on reversed sequence to get the beginning index
beginIndices= findMax(rev1[::-1], rev2[::-1])
end= endIndices[0][0],endIndices[0][1]
#refactor beginning indices to be correct in the context of the unaltered seq
begin= endIndices[0][0]-beginIndices[0][0], endIndices[0][1]-beginIndices[0][1] 
#print(begin, end)
org1=p[0]
org2=p[1]
#reduce strings to the found beginning and end points 
s1= org1[begin[0]:end[0]+1]
#print(s1)
s2= org2[begin[1]:end[1]+1]
#call local on the reduced sequences
l=locAL(s1,s2)
localQ1.write("high score of: "+ str(endIndices[1])+" "+"with length: "+ str(l[1])+"\n")
localQ1.close()
